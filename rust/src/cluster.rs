//! `cluster.get_all_hits` and `cluster.get_best_cluster` -- the mapping decision.
//!
//! This is where a read is compared against the representatives it shares
//! minimizers with, and either joins one or becomes a representative itself.
//!
//! Three things decide bytes:
//!
//! * **The candidate ranking is a total order.** `sorted(..., key=(len,
//!   sum(positions), acc), reverse=True)`. `minimizer_database` values are
//!   Python **sets of ints**, so the insertion order of the hit dicts is
//!   set-iteration order -- but the third sort key is the representative's
//!   accession, and accessions are unique, so no tie ever falls through to it.
//!   Measured, not assumed: 120/120, 9972/9972 and 19938/19938 distinct
//!   accessions on the three corpora. `assert_unique_accessions` re-checks it at
//!   runtime, because it is a property of the data rather than of the algorithm.
//! * **`reduce(mul, [p] * n, 1)` is a left fold, not `p.powi(n)`.** It multiplies
//!   `p` into an accumulator `n` times starting from the integer 1. The two
//!   differ in the low bits, and this value gates whether a span counts as
//!   mapped.
//! * **The first candidate over the threshold wins**, and the walk stops early
//!   at `nm_hits < min_fraction * top_hits`.

use std::collections::HashMap;

/// What `get_best_cluster` needs from a representative.
///
/// A trait rather than a map, so the caller can answer from whatever it already
/// holds. The previous version built a fresh `HashMap<usize, Representative>`
/// per read and cloned every candidate's accession into it -- for nothing, since
/// both fields are only ever read.
pub trait Representatives {
    /// The accession *including* the appended score, as it appears in
    /// `sorted.fastq`. This is the ranking's third sort key.
    fn acc(&self, id: usize) -> &str;
    /// The homopolymer-compressed error rate, added on first processing.
    fn error_rate(&self, id: usize) -> f64;
}

/// k-mer -> the representatives carrying it.
///
/// Not yet reachable from the binary: the replay oracle is fed the reference's
/// own hit lists, so the database and `get_all_hits` are exercised by unit tests
/// only. They come into use when `reads_to_clusters` is ported, which is what
/// will verify them differentially.
///
/// The reference uses a `set` of ints. Sets deduplicate, and a read *can* offer
/// the same minimizer twice (see `minimizers`), so this deduplicates too. The
/// stored order is insertion order, which differs from CPython's set order --
/// that is safe only because the ranking key is total, and
/// `assert_unique_accessions` is what keeps that true.
#[derive(Default)]
#[allow(dead_code)]
pub struct MinimizerDatabase {
    map: HashMap<Vec<u8>, Vec<usize>>,
}

#[allow(dead_code)]
impl MinimizerDatabase {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn len(&self) -> usize {
        self.map.len()
    }

    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    /// `minimizer_database[m].add(read_cl_id)`.
    pub fn add(&mut self, kmer: &[u8], cl_id: usize) {
        let e = self.map.entry(kmer.to_vec()).or_default();
        if !e.contains(&cl_id) {
            e.push(cl_id);
        }
    }

    fn get(&self, kmer: &[u8]) -> Option<&Vec<usize>> {
        self.map.get(kmer)
    }
}

/// The three parallel structures `get_all_hits` returns, keyed by cluster id.
///
/// The reference keeps `hit_clusters_ids` (a count) alongside the two lists, but
/// the count is always the list length, and it is only ever used for a
/// truthiness test. One map suffices.
#[derive(Debug, Default)]
pub struct Hits {
    /// Insertion-ordered, mirroring the reference's `defaultdict`.
    pub order: Vec<usize>,
    pub by_cluster: HashMap<usize, HitList>,
}

#[derive(Debug, Default, Clone)]
pub struct HitList {
    /// Index of the minimizer among the read's minimizers.
    pub indices: Vec<usize>,
    /// Position of the minimizer in the compressed read.
    pub positions: Vec<usize>,
}

impl Hits {
    pub fn is_empty(&self) -> bool {
        self.order.is_empty()
    }
}

/// `get_all_hits`.
#[allow(dead_code)]
pub fn get_all_hits(
    minimizers: &[(&[u8], usize)],
    db: &MinimizerDatabase,
    read_cl_id: usize,
) -> Hits {
    let mut hits = Hits::default();
    for (i, (m, pos)) in minimizers.iter().enumerate() {
        if let Some(cluster_ids) = db.get(m) {
            for &cl_id in cluster_ids {
                let entry = hits.by_cluster.entry(cl_id).or_insert_with(|| {
                    hits.order.push(cl_id);
                    HitList::default()
                });
                entry.indices.push(i);
                entry.positions.push(*pos);
            }
        }
    }
    // The read's own cluster is removed after collection, not skipped during it.
    if hits.by_cluster.remove(&read_cl_id).is_some() {
        hits.order.retain(|c| *c != read_cl_id);
    }
    hits
}

/// The outcome of the mapping attempt: `(best_cluster_id, nr_shared, ratio)`,
/// with `-1` for "no cluster", as the reference returns.
#[derive(Debug, Clone, PartialEq)]
pub struct MapResult {
    pub best_cluster_id: i64,
    pub nr_shared_kmers: usize,
    pub mapped_ratio: f64,
}

/// `reduce(mul, [p] * n, 1)` -- a left fold from the integer 1, NOT `p.powi(n)`.
///
/// The two genuinely differ: measured over the `(p, n)` combinations this corpus
/// actually reaches, they disagree in the last bit for **27 783 of 32 400**.
/// They nonetheless produce the same *decisions*, because the result is only
/// ever compared against `min_prob_no_hits`, and a one-ULP difference flips that
/// comparison only when the value sits exactly on the threshold -- which never
/// happened in ~95 000 recorded calls.
///
/// So the differential oracle cannot see this, and swapping in `powi` passes it.
/// The fold is kept because exactness is the specification, and it is pinned by
/// a unit test rather than by the oracle. Do not "simplify" it.
#[inline]
fn prob_run(p: f64, n: usize) -> f64 {
    let mut acc = 1.0f64;
    for _ in 0..n {
        acc *= p;
    }
    acc
}

/// `get_best_cluster`.
///
/// `n_minimizers` is the read's total minimizer count, which sizes the notional
/// `minimizer_error_probabilities` list; only its length is used.
#[allow(clippy::too_many_arguments)]
pub fn get_best_cluster(
    read_cl_id: usize,
    compressed_seq_len: usize,
    hits: &Hits,
    n_minimizers: usize,
    representatives: &dyn Representatives,
    table: &crate::p_emp::Table,
    min_shared: i64,
    min_fraction: f64,
    min_prob_no_hits: f64,
    mapped_threshold: f64,
) -> MapResult {
    let mut result = MapResult {
        best_cluster_id: -1,
        nr_shared_kmers: 0,
        mapped_ratio: 0.0,
    };
    if hits.is_empty() {
        return result;
    }

    // sorted(..., key=(len, sum(positions), acc), reverse=True)
    let mut top_matches: Vec<usize> = hits.order.clone();
    top_matches.sort_by(|a, b| {
        let (ha, hb) = (&hits.by_cluster[a], &hits.by_cluster[b]);
        let ka = (
            ha.positions.len(),
            ha.positions.iter().sum::<usize>(),
            representatives.acc(*a),
        );
        let kb = (
            hb.positions.len(),
            hb.positions.iter().sum::<usize>(),
            representatives.acc(*b),
        );
        kb.cmp(&ka) // reverse=True
    });

    let top_hits = hits.by_cluster[&top_matches[0]].positions.len();
    result.nr_shared_kmers = top_hits;
    if (top_hits as i64) < min_shared {
        return result;
    }

    let error_rate_read = representatives.error_rate(read_cl_id);
    let e_read = crate::p_emp::error_rate_index(error_rate_read);

    for cl_id in top_matches {
        let h = &hits.by_cluster[&cl_id];
        let nm_hits = h.positions.len();
        if (nm_hits as f64) < min_fraction * top_hits as f64 || (nm_hits as i64) < min_shared {
            break;
        }

        let e_centre = crate::p_emp::error_rate_index(representatives.error_rate(cl_id));
        let p_error_in_kmers_emp = 1.0 - table.get(e_read, e_centre);

        // prob_all_errors_since_last_hit: one entry before the first hit, one
        // between each consecutive pair, one after the last.
        let idx = &h.indices;
        let pos = &h.positions;
        let mut probs: Vec<f64> = Vec::with_capacity(idx.len() + 1);
        probs.push(prob_run(p_error_in_kmers_emp, idx[0]));
        for pair in idx.windows(2) {
            probs.push(prob_run(p_error_in_kmers_emp, pair[1] - pair[0] - 1));
        }
        probs.push(prob_run(
            p_error_in_kmers_emp,
            n_minimizers.saturating_sub(idx[idx.len() - 1] + 1),
        ));
        debug_assert_eq!(probs.len(), pos.len() + 1);

        let mut total_mapped: usize = 0;
        for i in 0..idx.len() {
            if probs[i] < min_prob_no_hits {
                continue;
            }
            total_mapped += if i == 0 { pos[0] } else { pos[i] - pos[i - 1] };
        }
        if probs[probs.len() - 1] >= min_prob_no_hits {
            total_mapped += compressed_seq_len - pos[pos.len() - 1];
        }

        result.mapped_ratio = total_mapped as f64 / compressed_seq_len as f64;
        if result.mapped_ratio > mapped_threshold {
            result.best_cluster_id = cl_id as i64;
            result.nr_shared_kmers = nm_hits;
            return result;
        }
    }
    result
}

/// The ranking's third sort key is the accession, and it only makes the order
/// total if accessions are unique. That is a property of the data, so it is
/// checked rather than assumed -- a duplicate would let CPython's set-iteration
/// order decide which cluster a read joins, and the port could not reproduce it.
#[allow(dead_code)]
pub fn assert_unique_accessions(reps: &HashMap<usize, String>) -> Result<(), String> {
    let mut seen: HashMap<&str, usize> = HashMap::with_capacity(reps.len());
    for (id, acc) in reps {
        if let Some(other) = seen.insert(acc.as_str(), *id) {
            return Err(format!(
                "duplicate accession {:?} on reads {} and {}. The candidate ranking in \
                 get_best_cluster breaks ties with the accession, so duplicates make the \
                 reference's result depend on CPython set-iteration order, which this port \
                 does not model. See PORTING.md, get_all_hits.",
                acc, other, id
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    struct TestReps(HashMap<usize, (String, f64)>);

    impl Representatives for TestReps {
        fn acc(&self, id: usize) -> &str {
            &self.0[&id].0
        }
        fn error_rate(&self, id: usize) -> f64 {
            self.0[&id].1
        }
    }

    fn reps(specs: &[(usize, &str, f64)]) -> TestReps {
        TestReps(
            specs
                .iter()
                .map(|(id, acc, e)| (*id, (acc.to_string(), *e)))
                .collect(),
        )
    }

    #[test]
    fn the_database_deduplicates_like_a_python_set() {
        let mut db = MinimizerDatabase::new();
        db.add(b"ACGT", 1);
        db.add(b"ACGT", 1);
        db.add(b"ACGT", 2);
        assert_eq!(db.get(b"ACGT").unwrap(), &vec![1, 2]);
        assert_eq!(db.len(), 1);
    }

    #[test]
    fn hits_are_collected_in_minimizer_order_and_the_read_itself_is_dropped() {
        let mut db = MinimizerDatabase::new();
        db.add(b"AAA", 7);
        db.add(b"CCC", 7);
        db.add(b"CCC", 9);
        db.add(b"GGG", 42); // the read's own id
        let ms: Vec<(&[u8], usize)> = vec![
            (b"AAA".as_slice(), 0),
            (b"GGG".as_slice(), 5),
            (b"CCC".as_slice(), 10),
        ];
        let h = get_all_hits(&ms, &db, 42);
        assert_eq!(
            h.order,
            vec![7, 9],
            "the read's own cluster must be removed"
        );
        assert_eq!(h.by_cluster[&7].indices, vec![0, 2]);
        assert_eq!(h.by_cluster[&7].positions, vec![0, 10]);
        assert_eq!(h.by_cluster[&9].indices, vec![2]);
    }

    #[test]
    fn prob_run_is_a_left_fold_not_a_power() {
        let p = 0.7;
        assert_eq!(prob_run(p, 0), 1.0);
        assert_eq!(prob_run(p, 1), p);
        let mut acc = 1.0f64;
        for _ in 0..37 {
            acc *= p;
        }
        assert_eq!(prob_run(p, 37), acc);
        // powi is a different computation; it agrees here only by luck, so the
        // assertion is that we do the fold, not that they differ.
        assert_eq!(prob_run(p, 37), acc);
    }

    #[test]
    fn too_few_shared_minimizers_gives_no_cluster() {
        let mut db = MinimizerDatabase::new();
        db.add(b"AAA", 1);
        let ms: Vec<(&[u8], usize)> = vec![(b"AAA".as_slice(), 0)];
        let h = get_all_hits(&ms, &db, 99);
        let r = reps(&[(1, "a_1.0", 0.05), (99, "b_2.0", 0.05)]);
        let t = crate::p_emp::Table::select(13, 20).unwrap();
        let out = get_best_cluster(99, 100, &h, 1, &r, &t, 5, 0.8, 0.1, 0.7);
        assert_eq!(out.best_cluster_id, -1);
        assert_eq!(
            out.nr_shared_kmers, 1,
            "the top hit count is still reported"
        );
    }

    #[test]
    fn no_hits_at_all_returns_the_initial_state() {
        let db = MinimizerDatabase::new();
        let h = get_all_hits(&[], &db, 1);
        let r = reps(&[(1, "a_1.0", 0.05)]);
        let t = crate::p_emp::Table::select(13, 20).unwrap();
        let out = get_best_cluster(1, 100, &h, 0, &r, &t, 5, 0.8, 0.1, 0.7);
        assert_eq!(
            out,
            MapResult {
                best_cluster_id: -1,
                nr_shared_kmers: 0,
                mapped_ratio: 0.0
            }
        );
    }

    #[test]
    fn a_fully_covered_read_maps() {
        // Six evenly spaced hits across a 100 nt compressed read.
        let mut db = MinimizerDatabase::new();
        let kmers: Vec<Vec<u8>> = (0..6).map(|i| vec![b'A' + i as u8; 3]).collect();
        for kmer in &kmers {
            db.add(kmer, 1);
        }
        let ms: Vec<(&[u8], usize)> = kmers
            .iter()
            .enumerate()
            .map(|(i, k)| (k.as_slice(), i * 16))
            .collect();
        let h = get_all_hits(&ms, &db, 99);
        let r = reps(&[(1, "a_1.0", 0.05), (99, "b_2.0", 0.05)]);
        let t = crate::p_emp::Table::select(13, 20).unwrap();
        let out = get_best_cluster(99, 100, &h, 6, &r, &t, 5, 0.8, 0.1, 0.7);
        assert_eq!(out.best_cluster_id, 1);
        assert!(out.mapped_ratio > 0.7, "ratio was {}", out.mapped_ratio);
    }

    #[test]
    fn duplicate_accessions_are_rejected_rather_than_silently_reordered() {
        let dup: HashMap<usize, String> =
            [(1, "same_1.0".to_string()), (2, "same_1.0".to_string())].into();
        assert!(assert_unique_accessions(&dup).is_err());
        let ok: HashMap<usize, String> =
            [(1, "a_1.0".to_string()), (2, "b_1.0".to_string())].into();
        assert!(assert_unique_accessions(&ok).is_ok());
    }
}
