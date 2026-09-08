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
use rustc_hash::FxHashMap;

/// A read as the sweep sees it: the reference's
/// `(read_cl_id, prev_batch_index, acc, seq, qual, score)`.
#[derive(Clone)]
pub struct SweepRead {
    pub id: usize,
    pub prev_batch_index: i64,
    /// Shared, like `qual` below: the accession is held by this read,
    /// by its `ReadInfo`, and again in the cluster's member list, so cloning the
    /// `String` meant three resident copies of every accession.
    pub acc: std::sync::Arc<str>,
    /// Shared with the `ReadInfo` for the same read rather than cloned into it.
    /// Every read starts as its own representative, so cloning meant a second
    /// resident copy of every base and quality score -- 129 MB of droso_100k's
    /// heap peak. `Arc` keeps the existing ownership structure, which matters
    /// because `ReadInfo` is moved between passes in parallel mode.
    /// 2-bit packed; see `packed`. Shared with this read's `ReadInfo`.
    pub seq: crate::packed::PackedSeq,
    /// `compressed_error_rate(seq, qual)`, computed at load.
    ///
    /// The quality string itself is **not** kept: it was 716 MB of a 1419 MB
    /// live heap on SIRV_real_full, half of everything, and `cluster.rs` never
    /// reads it. Everything the clustering needs from it is these two numbers,
    /// and the only consumer of the string is the origins writer, which streams
    /// `sorted.fastq` again for the surviving representatives.
    ///
    /// Precomputed but revealed lazily: `reads_to_clusters` copies this into
    /// `ReadInfo::error_rate` at exactly the point the reference computes it, so
    /// a read that never reaches that step still reports `nan`. Setting it
    /// eagerly would change those reads' output.
    pub hp_error_rate: Option<f64>,
    /// `expected_errors(qual) / seq.len() as f64`, computed at load. Written as
    /// that exact expression so the f64 is bit-identical to the reference's.
    pub err_per_base: f64,
    pub score: f64,
}

/// A representative. `error_rate` is `None` until step 2 fills it in, mirroring
/// the reference's 6-tuple that becomes a 7-tuple.
#[derive(Clone)]
pub struct ReadInfo {
    pub id: usize,
    pub batch_index: i64,
    /// Shared with the read's `SweepRead`; see the note there.
    pub acc: std::sync::Arc<str>,
    /// Shared with the read's `SweepRead`; see the note there.
    pub seq: crate::packed::PackedSeq,
    /// See `SweepRead::err_per_base`.
    pub err_per_base: f64,
    pub score: f64,
    pub error_rate: Option<f64>,
}

/// An insertion-ordered map, because the reference's dicts are.
///
/// `order` is load-bearing, but not in `main`: that sort is total -- `(size,
/// score, cluster id)` -- so insertion order cannot reach the final output.
/// The consumer is `parallelize::render_intermediate`, which sorts by cluster
/// size *only* with a stable `sort_by`, so insertion order breaks its ties and
/// reaches `pre_clusters.csv` and `cluster_origins.csv`.
///
/// Members are **read ids**, not accessions. The accession is recovered by
/// indexing the sorted read array, which is how the reference's own output loop
/// finds it too; holding one `Arc<str>` per member cost 20 MB of pointers on
/// SIRV_real_full where 5 MB of `u32` does.
#[derive(Default, Clone)]
pub struct OrderedClusters {
    pub order: Vec<usize>,
    pub map: FxHashMap<usize, Vec<u32>>,
}

impl OrderedClusters {
    pub fn insert(&mut self, id: usize, members: Vec<u32>) {
        if self.map.insert(id, members).is_none() {
            self.order.push(id);
        }
    }
    pub fn remove(&mut self, id: usize) -> Option<Vec<u32>> {
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
    pub fn iter(&self) -> impl Iterator<Item = (usize, &Vec<u32>)> {
        self.order.iter().map(move |i| (*i, &self.map[i]))
    }
}

/// Answers `get_best_cluster_block_align`'s questions by reference.
///
/// The previous version cloned each candidate's full sequence *and* quality
/// string per comparison -- roughly 12 KB per candidate on Drosophila reads, for
/// data that is only read.
struct RepSeqs<'a>(&'a FxHashMap<usize, ReadInfo>);

impl blockalign::AlignSource for RepSeqs<'_> {
    fn seq_into(&self, id: usize, out: &mut Vec<u8>) {
        self.0[&id].seq.unpack_into(out);
    }
    fn seq_len(&self, id: usize) -> usize {
        self.0[&id].seq.len()
    }
    fn err_per_base(&self, id: usize) -> f64 {
        self.0[&id].err_per_base
    }
    fn acc(&self, id: usize) -> &str {
        &self.0[&id].acc
    }
}

/// Answers `get_best_cluster`'s questions straight out of the sweep's own map.
///
/// This replaced a per-read `HashMap<usize, Representative>` that cloned every
/// candidate's accession. Behaviour-neutral: both fields are read-only.
struct RepMap<'a>(&'a FxHashMap<usize, ReadInfo>);

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

/// Per-stage wall clock, for `ISONCLUST_PROFILE=1`.
///
/// Explicit timers rather than a sampling profiler: at `--release` the stage
/// functions inline into `reads_to_clusters`, so `sample` attributes 94-99% of
/// everything to one symbol and cannot separate them. Method point 5 also says
/// to *remove* sub-stage instrumentation after reading it -- these are cheap
/// (one `Instant::now` per read per stage, not per inner loop) and off unless
/// asked for, but they should not outlive their usefulness.
#[derive(Default, Clone, Copy)]
pub struct StageTimes {
    pub minimizers: std::time::Duration,
    pub error_rate: std::time::Duration,
    pub hits: std::time::Duration,
    pub mapping: std::time::Duration,
    pub alignment: std::time::Duration,
    pub db_insert: std::time::Duration,
}

impl StageTimes {
    pub fn add(&mut self, o: &StageTimes) {
        self.minimizers += o.minimizers;
        self.error_rate += o.error_rate;
        self.hits += o.hits;
        self.mapping += o.mapping;
        self.alignment += o.alignment;
        self.db_insert += o.db_insert;
    }
    pub fn report(&self, label: &str) {
        let total = self.minimizers
            + self.error_rate
            + self.hits
            + self.mapping
            + self.alignment
            + self.db_insert;
        let t = total.as_secs_f64().max(1e-9);
        eprintln!(
            "  stage profile ({label}), {:.2}s accounted for:",
            total.as_secs_f64()
        );
        for (name, d) in [
            ("alignment (parasail)", self.alignment),
            ("mapping decision", self.mapping),
            ("minimizers", self.minimizers),
            ("hit collection", self.hits),
            ("compressed error rate", self.error_rate),
            ("database insert", self.db_insert),
        ] {
            eprintln!(
                "    {:<24} {:7.2}s  {:5.1}%",
                name,
                d.as_secs_f64(),
                100.0 * d.as_secs_f64() / t
            );
        }
    }
}

/// What one sweep returns, mirroring the reference's
/// `{new_batch_index: (clusters, representatives, minimizer_database, new_batch_index)}`.
pub struct SweepResult {
    pub clusters: OrderedClusters,
    pub representatives: FxHashMap<usize, ReadInfo>,
    pub db: MinimizerDatabase,
    pub batch_index: i64,
    pub mapped_passed: usize,
    pub aln_passed: usize,
    pub aln_called: usize,
    pub skipped_short: usize,
    pub times: StageTimes,
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
    mut reps: FxHashMap<usize, ReadInfo>,
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

    // Which reads have been merged away, kept only for the one-level-deep
    // invariant the reference relies on. Debug builds check it; release builds
    // still need the set, so it stays.
    let mut merged: std::collections::HashSet<usize> = std::collections::HashSet::new();
    // One unpacking buffer for the whole sweep; see `packed`.
    let mut seq_buf: Vec<u8> = Vec::new();
    // One hit collector for the whole sweep, reset per read; see `cluster::Hits`.
    let mut hits_buf = cluster::Hits::default();
    let mut out_mapped = 0usize;
    let mut out_aln_passed = 0usize;
    let mut out_aln_called = 0usize;
    let mut skipped_short = 0usize;
    let mut times = StageTimes::default();
    let profiling = std::env::var("ISONCLUST_PROFILE").is_ok();
    macro_rules! timed {
        ($field:ident, $body:expr) => {{
            if profiling {
                let t = std::time::Instant::now();
                let v = $body;
                times.$field += t.elapsed();
                v
            } else {
                $body
            }
        }};
    }

    for r in sorted_reads {
        let read_cl_id = r.id;

        // Every read starts as its own representative, as the reference's
        // initialisation does -- but created here, on first sight, rather than
        // pre-populated for the whole corpus before the sweep. A read that maps
        // is removed again immediately below, so the two maps hold roughly the
        // surviving clusters instead of every read: on SIRV_real_full that is
        // 579 entries rather than 1 295 814, which is 260 MB of table and
        // per-read Vec that used to be built and then thrown away.
        //
        // It has to happen before the early exits, not after: a read skipped for
        // being too short still reaches the output as its own cluster, so it
        // needs its entry even though it is never processed.
        //
        // `or_insert_with` also means a caller that pre-populates (parallelize,
        // and later passes carrying representatives forward) is unaffected, and
        // that a carried-over `error_rate` is never overwritten.
        reps.entry(read_cl_id).or_insert_with(|| ReadInfo {
            id: r.id,
            batch_index: r.prev_batch_index,
            acc: r.acc.clone(),
            seq: r.seq.clone(),
            err_per_base: r.err_per_base,
            score: r.score,
            error_rate: None,
        });
        if !clusters.map.contains_key(&read_cl_id) {
            clusters.insert(read_cl_id, vec![read_cl_id as u32]);
        }

        if r.prev_batch_index == lowest_batch_index {
            if let Some(info) = reps.get_mut(&read_cl_id) {
                info.batch_index = new_batch_index;
            }
            continue;
        }

        // 1. compress and take minimizers
        //
        // Unpacked into a buffer reused across every read, so packing the
        // sequences costs no per-read allocation here.
        r.seq.unpack_into(&mut seq_buf);
        let hpol = crate::sorting::homopolymer_compress(&seq_buf);
        if hpol.len() < p.k {
            skipped_short += 1;
            continue;
        }
        let ms = timed!(minimizers, minimizers::get_kmer_minimizers(&hpol, p.k, p.w));

        // 2. the compressed error rate, unless a previous pass already did it
        {
            let info = reps
                .get_mut(&read_cl_id)
                .expect("read has a representative");
            if info.error_rate.is_some() {
                info.batch_index = new_batch_index;
            } else {
                info.batch_index = new_batch_index;
                // Computed at load; see `SweepRead::hp_error_rate`.
                info.error_rate = r.hp_error_rate;
            }
        }

        // 3. hits
        timed!(
            hits,
            cluster::get_all_hits(&ms, &db, read_cl_id, &mut hits_buf)
        );

        // 4. map
        let m = timed!(
            mapping,
            cluster::get_best_cluster(
                read_cl_id,
                hpol.len(),
                &hits_buf,
                ms.len(),
                &RepMap(&reps),
                table,
                p.min_shared,
                p.min_fraction,
                p.min_prob_no_hits,
                p.mapped_threshold,
            )
        );
        if m.best_cluster_id >= 0 {
            out_mapped += 1;
        }

        // 5. align
        let a_id = if m.best_cluster_id < 0 && (m.nr_shared_kmers as i64) >= p.min_shared {
            out_aln_called += 1;
            let a = timed!(
                alignment,
                blockalign::get_best_cluster_block_align(
                    read_cl_id,
                    &hits_buf,
                    &RepSeqs(&reps),
                    p.k,
                    p.aligned_threshold,
                )
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
            // Applied now rather than collected and replayed after the loop.
            // Safe because nothing inside this loop reads cluster membership or
            // size -- only `hpol.len()` and `ms.len()` -- and because a merge
            // target is always a representative while representatives are never
            // merge sources, so a read removed here can never be needed again:
            // candidates come from the minimizer database, which only ever holds
            // representatives. The order targets receive members in is the loop
            // order either way, so the member lists are unchanged.
            let target = best as usize;
            debug_assert!(
                !merged.contains(&target),
                "merge target {target} is itself merged; the reference assumes this cannot happen"
            );
            merged.insert(read_cl_id);
            let moved = clusters
                .remove(read_cl_id)
                .expect("a source is merged once");
            clusters
                .map
                .get_mut(&target)
                .expect("a target is never merged away")
                .extend(moved);
            reps.remove(&read_cl_id);
        } else {
            timed!(db_insert, {
                for (mn, _) in &ms {
                    db.add(mn, read_cl_id);
                }
            });
        }
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
        times,
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
            c.insert(i, vec![i as u32]);
        }
        c.remove(2);
        assert_eq!(c.order, vec![0, 1, 3, 4]);
        c.insert(9, vec![9u32]);
        assert_eq!(c.order, vec![0, 1, 3, 4, 9]);
    }
}
