//! `cluster.reads_to_clusters` -- the greedy sweep.
//!
//! Reads are visited highest-score first. Each is compared against the
//! representatives it shares minimizers with; if it joins one it contributes
//! nothing further, and if it does not it becomes a representative and **all of
//! its minimizers enter the database**. The database only grows, so a read is
//! only ever compared against representatives that already exist -- which is why
//! the score order has to be exact.
//!
//! The signature carries the multiprocessing machinery even in single-core
//! mode, because the reference's does: `parallel_clustering` calls this same
//! function with a pre-populated cluster map, a carried-over minimizer database,
//! and a batch index, and the sweep's first act is to **skip every read whose
//! previous batch index equals the lowest in the batch** -- those are the reads
//! that built the database being reused.
//!
//! Order of business per read, matching the reference's numbered comments:
//!
//! 1. homopolymer-compress and take minimizers (short reads are skipped)
//! 2. compute the compressed-read error rate, unless it is already there --
//!    the reference branches on `len(representatives[id]) == 7`
//! 3. collect hits against the database
//! 4. try to map
//! 5. if mapping failed but at least `--min_shared` minimizers were shared, align
//! 6. record the assignment, or become a representative and add the minimizers
//! 7. reassign: move every recorded read into its target's cluster

use crate::blockalign;
use crate::cluster::{self, MinimizerDatabase, Representatives};
use crate::minimizers;
use std::collections::HashMap;

/// A read as the sweep sees it: the reference's
/// `(read_cl_id, prev_batch_index, acc, seq, qual, score)`.
#[derive(Clone)]
pub struct SweepRead {
    pub id: usize,
    pub prev_batch_index: i64,
    pub acc: String,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub score: f64,
}

/// A representative. `error_rate` is `None` until step 2 fills it in, mirroring
/// the reference's 6-tuple that becomes a 7-tuple.
#[derive(Clone)]
pub struct ReadInfo {
    pub id: usize,
    pub batch_index: i64,
    pub acc: String,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub score: f64,
    pub error_rate: Option<f64>,
}

/// An insertion-ordered map, because the reference's dicts are and the order
/// reaches the output through the stable sort in `main`.
#[derive(Default, Clone)]
pub struct OrderedClusters {
    pub order: Vec<usize>,
    pub map: HashMap<usize, Vec<String>>,
}

impl OrderedClusters {
    pub fn insert(&mut self, id: usize, accs: Vec<String>) {
        if self.map.insert(id, accs).is_none() {
            self.order.push(id);
        }
    }
    pub fn remove(&mut self, id: usize) -> Option<Vec<String>> {
        let v = self.map.remove(&id);
        if v.is_some() {
            self.order.retain(|x| *x != id);
        }
        v
    }
    /// Used by the tests and by callers that report cluster counts.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.order.len()
    }
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.order.is_empty()
    }
    pub fn iter(&self) -> impl Iterator<Item = (usize, &Vec<String>)> {
        self.order.iter().map(move |i| (*i, &self.map[i]))
    }
}

/// Answers `get_best_cluster_block_align`'s questions by reference.
///
/// The previous version cloned each candidate's full sequence *and* quality
/// string per comparison -- roughly 12 KB per candidate on Drosophila reads, for
/// data that is only read.
struct RepSeqs<'a>(&'a HashMap<usize, ReadInfo>);

impl blockalign::AlignSource for RepSeqs<'_> {
    fn seq_qual(&self, id: usize) -> (&[u8], &[u8]) {
        let x = &self.0[&id];
        (&x.seq, &x.qual)
    }
    fn acc(&self, id: usize) -> &str {
        &self.0[&id].acc
    }
}

/// Answers `get_best_cluster`'s questions straight out of the sweep's own map.
///
/// This replaced a per-read `HashMap<usize, Representative>` that cloned every
/// candidate's accession. Behaviour-neutral: both fields are read-only.
struct RepMap<'a>(&'a HashMap<usize, ReadInfo>);

impl Representatives for RepMap<'_> {
    fn acc(&self, id: usize) -> &str {
        &self.0[&id].acc
    }
    fn error_rate(&self, id: usize) -> f64 {
        self.0[&id].error_rate.unwrap_or(f64::NAN)
    }
}

/// Tunables, passed straight through from the CLI.
#[derive(Clone, Copy)]
pub struct SweepParams {
    pub k: usize,
    pub w: usize,
    pub min_shared: i64,
    pub min_fraction: f64,
    pub min_prob_no_hits: f64,
    pub mapped_threshold: f64,
    pub aligned_threshold: f64,
}

/// What one sweep returns, mirroring the reference's
/// `{new_batch_index: (clusters, representatives, minimizer_database, new_batch_index)}`.
pub struct SweepResult {
    pub clusters: OrderedClusters,
    pub representatives: HashMap<usize, ReadInfo>,
    pub db: MinimizerDatabase,
    pub batch_index: i64,
    pub mapped_passed: usize,
    pub aln_passed: usize,
    pub aln_called: usize,
    pub skipped_short: usize,
}

/// The compressed quality string: one character per homopolymer run, the best of
/// the run.
///
/// `min(qual[start:start+len], key=phred)` picks the *lowest error probability*,
/// i.e. the highest quality, and Python's `min` returns the first such character
/// on a tie.
pub fn compressed_quality(seq: &[u8], qual: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    let mut i = 0usize;
    while i < seq.len() {
        let mut j = i + 1;
        while j < seq.len() && seq[j] == seq[i] {
            j += 1;
        }
        let run = &qual[i.min(qual.len())..j.min(qual.len())];
        if let Some(best) = run.iter().copied().reduce(|a, b| {
            if crate::phred::capped(b) < crate::phred::capped(a) {
                b
            } else {
                a
            }
        }) {
            out.push(best);
        }
        i = j;
    }
    out
}

/// The homopolymer-compressed error rate step 2 appends.
pub fn compressed_error_rate(seq: &[u8], qual: &[u8]) -> Option<f64> {
    let qc = compressed_quality(seq, qual);
    if qc.is_empty() {
        return None;
    }
    Some(blockalign::expected_errors(&qc) / qc.len() as f64)
}

/// `reads_to_clusters(clusters, representatives, sorted_reads, p_emp_probs,
/// minimizer_database, new_batch_index, args)`.
pub fn reads_to_clusters(
    mut clusters: OrderedClusters,
    mut reps: HashMap<usize, ReadInfo>,
    sorted_reads: &[SweepRead],
    mut db: MinimizerDatabase,
    new_batch_index: i64,
    table: &crate::p_emp::Table,
    p: SweepParams,
) -> SweepResult {
    // The reads that built the database being reused are skipped rather than
    // re-clustered. On the first pass every prev index is 0, so `max(1, min)` is
    // 1 and nothing matches -- "Saved: 0 iterations."
    let lowest_batch_index = sorted_reads
        .iter()
        .map(|r| r.prev_batch_index)
        .min()
        .unwrap_or(0)
        .max(1);

    let mut assignment: Vec<(usize, usize)> = Vec::new();
    let mut out_mapped = 0usize;
    let mut out_aln_passed = 0usize;
    let mut out_aln_called = 0usize;
    let mut skipped_short = 0usize;

    for r in sorted_reads {
        let read_cl_id = r.id;

        if r.prev_batch_index == lowest_batch_index {
            if let Some(info) = reps.get_mut(&read_cl_id) {
                info.batch_index = new_batch_index;
            }
            continue;
        }

        // 1. compress and take minimizers
        let hpol = crate::sorting::homopolymer_compress(&r.seq);
        if hpol.len() < p.k {
            skipped_short += 1;
            continue;
        }
        let ms = minimizers::get_kmer_minimizers(&hpol, p.k, p.w);

        // 2. the compressed error rate, unless a previous pass already did it
        {
            let info = reps
                .get_mut(&read_cl_id)
                .expect("read has a representative");
            if info.error_rate.is_some() {
                info.batch_index = new_batch_index;
            } else {
                info.batch_index = new_batch_index;
                info.error_rate = compressed_error_rate(&r.seq, &r.qual);
            }
        }

        // 3. hits
        let hits = cluster::get_all_hits(&ms, &db, read_cl_id);

        // 4. map
        let m = cluster::get_best_cluster(
            read_cl_id,
            hpol.len(),
            &hits,
            ms.len(),
            &RepMap(&reps),
            table,
            p.min_shared,
            p.min_fraction,
            p.min_prob_no_hits,
            p.mapped_threshold,
        );
        if m.best_cluster_id >= 0 {
            out_mapped += 1;
        }

        // 5. align
        let a_id = if m.best_cluster_id < 0 && (m.nr_shared_kmers as i64) >= p.min_shared {
            out_aln_called += 1;
            let a = blockalign::get_best_cluster_block_align(
                read_cl_id,
                &hits,
                &RepSeqs(&reps),
                p.k,
                p.aligned_threshold,
            );
            if a.best_cluster_id >= 0 {
                out_aln_passed += 1;
            }
            a.best_cluster_id
        } else {
            -1
        };

        // 6. assign, or become a representative
        let best = m.best_cluster_id.max(a_id);
        if best >= 0 {
            assignment.push((read_cl_id, best as usize));
        } else {
            for (mn, _) in &ms {
                db.add(mn, read_cl_id);
            }
        }
    }

    // 7. reassign. One level deep: a merge target is always a representative and
    //    representatives are never merge sources.
    let sources: std::collections::HashSet<usize> = assignment.iter().map(|(r, _)| *r).collect();
    for (_, target) in &assignment {
        debug_assert!(
            !sources.contains(target),
            "merge target {target} is itself merged; the reference assumes this cannot happen"
        );
    }
    for (read_cl_id, target) in assignment {
        let moved = clusters
            .remove(read_cl_id)
            .expect("a source is merged once");
        clusters
            .map
            .get_mut(&target)
            .expect("a target is never merged away")
            .extend(moved);
        reps.remove(&read_cl_id);
    }

    SweepResult {
        clusters,
        representatives: reps,
        db,
        batch_index: new_batch_index,
        mapped_passed: out_mapped,
        aln_passed: out_aln_passed,
        aln_called: out_aln_called,
        skipped_short,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compressed_quality_takes_the_best_of_each_run() {
        assert_eq!(compressed_quality(b"AACCG", b"I!#5J"), b"I5J".to_vec());
    }

    #[test]
    fn compressed_quality_breaks_ties_on_the_first_character() {
        assert_eq!(compressed_quality(b"AA", b"II"), b"I".to_vec());
    }

    #[test]
    fn no_homopolymers_leaves_the_quality_alone() {
        assert_eq!(compressed_quality(b"ACGT", b"IJKL"), b"IJKL".to_vec());
    }

    #[test]
    fn compressed_error_rate_uses_the_capped_table() {
        assert_eq!(
            compressed_error_rate(b"A", b"!").expect("non-empty"),
            0.79433
        );
    }

    #[test]
    fn ordered_clusters_keeps_insertion_order_through_removals() {
        let mut c = OrderedClusters::default();
        for i in 0..5 {
            c.insert(i, vec![format!("r{i}")]);
        }
        c.remove(2);
        assert_eq!(c.order, vec![0, 1, 3, 4]);
        c.insert(9, vec!["r9".into()]);
        assert_eq!(c.order, vec![0, 1, 3, 4, 9]);
    }
}
