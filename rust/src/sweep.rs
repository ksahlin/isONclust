//! `cluster.reads_to_clusters` -- the greedy sweep.
//!
//! Reads are visited highest-score first. Each is compared against the
//! representatives it shares minimizers with; if it joins one it contributes
//! nothing further, and if it does not it becomes a representative and **all of
//! its minimizers enter the database**. The database only grows, so a read is
//! only ever compared against representatives that already exist -- which is why
//! the score order has to be exact, and why this cannot be parallelised without
//! changing the answer.
//!
//! Order of business per read, matching the reference's numbered comments:
//!
//! 1. homopolymer-compress and take minimizers (short reads are skipped here,
//!    with a message on stdout)
//! 2. compute the compressed-read error rate, appending it to the
//!    representative tuple -- the reference branches on `len(...) == 7` to
//!    decide whether this is already done
//! 3. collect hits against the database
//! 4. try to map
//! 5. if mapping failed but at least `--min_shared` minimizers were shared, try
//!    to align
//! 6. record the assignment, or become a representative and add the minimizers
//! 7. reassign: move every recorded read into its target's cluster
//!
//! Step 7 is one level deep. Merge targets are always representatives, and
//! representatives never appear as merge sources, so there are no chains to
//! follow -- asserted rather than assumed, because a chain would silently
//! corrupt the output.

use crate::blockalign;
use crate::cluster::{self, MinimizerDatabase, Representative};
use crate::minimizers;
use crate::sorting::Scored;
use std::collections::HashMap;

/// What the sweep produces: the clusters, and the representative of each.
pub struct SweepResult {
    /// cluster id -> the accessions in it, in the order they were added.
    pub clusters: Vec<(usize, Vec<String>)>,
    pub representatives: HashMap<usize, ReadInfo>,
    pub mapped_passed: usize,
    pub aln_passed: usize,
    pub aln_called: usize,
    pub skipped_short: Vec<usize>,
}

#[derive(Clone)]
pub struct ReadInfo {
    pub acc: String,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub score: f64,
    pub error_rate: f64,
}

/// The compressed quality string: one character per homopolymer run, the best
/// of the run.
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

/// The homopolymer-compressed error rate that step 2 appends.
pub fn compressed_error_rate(seq: &[u8], qual: &[u8]) -> Option<f64> {
    let qc = compressed_quality(seq, qual);
    if qc.is_empty() {
        return None;
    }
    let poisson_mean = blockalign::expected_errors(&qc);
    Some(poisson_mean / qc.len() as f64)
}

/// `reads_to_clusters` for the single-core path (`--t 1`).
#[allow(clippy::too_many_arguments)]
pub fn reads_to_clusters(
    reads: &[Scored],
    k: usize,
    w: usize,
    table: &crate::p_emp::Table,
    min_shared: i64,
    min_fraction: f64,
    min_prob_no_hits: f64,
    mapped_threshold: f64,
    aligned_threshold: f64,
) -> SweepResult {
    let n = reads.len();
    // clusters[i] starts as [reads[i].acc]; insertion order is 0..n, and the
    // output ordering depends on it, so a Vec of Options preserves it exactly
    // while allowing the deletions step 7 makes.
    let mut clusters: Vec<Option<Vec<String>>> =
        reads.iter().map(|r| Some(vec![r.acc.clone()])).collect();
    let mut reps: HashMap<usize, ReadInfo> = HashMap::with_capacity(n);
    for (i, r) in reads.iter().enumerate() {
        reps.insert(
            i,
            ReadInfo {
                acc: r.acc.clone(),
                seq: r.seq.as_bytes().to_vec(),
                qual: r.qual.as_bytes().to_vec(),
                score: r.score,
                error_rate: f64::NAN, // filled in at step 2
            },
        );
    }

    let mut db = MinimizerDatabase::new();
    let mut assignment: Vec<(usize, usize)> = Vec::new(); // (read, target), in order
    let mut out = SweepResult {
        clusters: Vec::new(),
        representatives: HashMap::new(),
        mapped_passed: 0,
        aln_passed: 0,
        aln_called: 0,
        skipped_short: Vec::new(),
    };

    for read_cl_id in 0..n {
        let seq = reps[&read_cl_id].seq.clone();
        let qual = reps[&read_cl_id].qual.clone();

        // 1. compress and take minimizers
        let hpol = crate::sorting::homopolymer_compress(&seq);
        if hpol.len() < k {
            out.skipped_short.push(read_cl_id);
            continue;
        }
        let ms = minimizers::get_kmer_minimizers(&hpol, k, w);

        // 2. the compressed error rate
        let er = compressed_error_rate(&seq, &qual).unwrap_or(f64::NAN);
        reps.get_mut(&read_cl_id).expect("present").error_rate = er;

        // 3. hits
        let hits = cluster::get_all_hits(&ms, &db, read_cl_id);

        // 4. map
        let rep_view: HashMap<usize, Representative> = hits
            .order
            .iter()
            .chain(std::iter::once(&read_cl_id))
            .map(|id| {
                let r = &reps[id];
                (
                    *id,
                    Representative {
                        acc: r.acc.clone(),
                        error_rate: r.error_rate,
                    },
                )
            })
            .collect();
        let m = cluster::get_best_cluster(
            read_cl_id,
            hpol.len(),
            &hits,
            ms.len(),
            &rep_view,
            table,
            min_shared,
            min_fraction,
            min_prob_no_hits,
            mapped_threshold,
        );
        if m.best_cluster_id >= 0 {
            out.mapped_passed += 1;
        }

        // 5. align, only when mapping failed but enough minimizers were shared
        let a_id = if m.best_cluster_id < 0 && (m.nr_shared_kmers as i64) >= min_shared {
            out.aln_called += 1;
            let seqs = |id: usize| -> (Vec<u8>, Vec<u8>) {
                let r = &reps[&id];
                (r.seq.clone(), r.qual.clone())
            };
            let accs = |id: usize| -> String { reps[&id].acc.clone() };
            let a = blockalign::get_best_cluster_block_align(
                read_cl_id,
                &hits,
                &seqs,
                &accs,
                k,
                aligned_threshold,
            );
            if a.best_cluster_id >= 0 {
                out.aln_passed += 1;
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

    // 7. reassign. One level deep, asserted: a merge target is always a
    //    representative, and representatives are never merge sources.
    let sources: std::collections::HashSet<usize> = assignment.iter().map(|(r, _)| *r).collect();
    for (_, target) in &assignment {
        debug_assert!(
            !sources.contains(target),
            "merge target {target} is itself merged; the reference assumes this cannot happen"
        );
    }
    for (read_cl_id, target) in assignment {
        let moved = clusters[read_cl_id]
            .take()
            .expect("a source is merged once");
        clusters[target]
            .as_mut()
            .expect("a target is never merged away")
            .extend(moved);
        reps.remove(&read_cl_id);
    }

    out.clusters = clusters
        .into_iter()
        .enumerate()
        .filter_map(|(i, c)| c.map(|v| (i, v)))
        .collect();
    out.representatives = reps;
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compressed_quality_takes_the_best_of_each_run() {
        // runs: AA CC G  -> best of "I!" = 'I', best of "#5" = '5', "J" = 'J'
        let q = compressed_quality(b"AACCG", b"I!#5J");
        assert_eq!(q, b"I5J".to_vec());
    }

    #[test]
    fn compressed_quality_breaks_ties_on_the_first_character() {
        // both 'I': min returns the first
        assert_eq!(compressed_quality(b"AA", b"II"), b"I".to_vec());
    }

    #[test]
    fn no_homopolymers_leaves_the_quality_alone() {
        assert_eq!(compressed_quality(b"ACGT", b"IJKL"), b"IJKL".to_vec());
    }

    #[test]
    fn compressed_error_rate_uses_the_capped_table() {
        // '!' is phred 0, capped from 1.0 to 0.79433
        let e = compressed_error_rate(b"A", b"!").expect("non-empty");
        assert_eq!(e, 0.79433);
    }
}
