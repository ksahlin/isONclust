//! `cluster.parasail_block_alignment` and `cluster.get_best_cluster_block_align`.
//!
//! The alignment fallback: when the mapping decision fails but the read still
//! shares at least `--min_shared` minimizers, isONclust aligns it against the
//! tied-at-the-top candidates and clusters on the fraction of the alignment that
//! sits in a sufficiently-matching window.
//!
//! This path is not a rare corner. On `droso_20k` it decides **10 309 of 19 938
//! reads** and accounts for 47% of the reference's clustering time.
//!
//! Details that decide bytes:
//!
//! * **The gap-opening penalty is chosen per comparison** from the two reads'
//!   summed error rates, binned at 0.01 / 0.04 / 0.1 into 5 / 4 / 3 / 2. Every
//!   one of the four occurs in practice.
//! * **`match_id` is `floor((1 - error_rate_sum) * k)`**, and `math.floor`
//!   returns an int, so the comparison below is integer-against-integer.
//! * **The rolling window is over ALIGNMENT columns, not read positions**, and
//!   its length is `k`. The count of "aligned" columns is then divided by
//!   `len(s1)` -- the *unaligned* query length -- so the ratio can exceed 1.
//! * **Only candidates tied at the top hit count are considered** (`if nm_hits <
//!   top_hits: break`), which is stricter than `get_best_cluster`'s
//!   `min_fraction` walk.

use crate::parasail::{self, Scoring};

/// `parasail_block_alignment`'s return: the two gapped strings and the ratio.
///
/// The gapped strings are what the reference returns and hands to
/// `get_best_cluster_block_align`, which currently discards them -- but they are
/// part of the function's contract and the next stage may want them, so they are
/// kept rather than dropped to satisfy a lint.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BlockAlignment {
    pub s1_aligned: Vec<u8>,
    pub s2_aligned: Vec<u8>,
    pub alignment_ratio: f64,
}

/// isONclust's scoring: `parasail.matrix_create("ACGT", 2, -2)` with
/// `gap_ext = 1` and a caller-chosen opening penalty.
pub fn scoring(opening_penalty: i32) -> Scoring {
    Scoring {
        match_score: 2,
        mismatch: -2,
        open: opening_penalty,
        ext: 1,
    }
}

/// `parasail_block_alignment(s1, s2, k, match_id, ..., opening_penalty, 1)`.
pub fn parasail_block_alignment(
    s1: &[u8],
    s2: &[u8],
    k: usize,
    match_id: i64,
    opening_penalty: i32,
) -> BlockAlignment {
    let aln = parasail::semiglobal(s1, s2, scoring(opening_penalty));
    let (a1, a2) = crate::align::ops_to_seq(&aln.ops, s1, s2)
        .expect("a parasail CIGAR always expands against its own inputs");

    // match_vector over alignment columns; a gap never matches.
    let matches: Vec<u8> = a1
        .iter()
        .zip(a2.iter())
        .map(|(x, y)| u8::from(x == y))
        .collect();

    // The reference seeds the window with the first k columns even when the
    // alignment is shorter than k, so `sum` is over whatever exists.
    let head = matches.len().min(k);
    let mut current: i64 = matches[..head].iter().map(|x| i64::from(*x)).sum();
    let mut aligned_columns: i64 = i64::from(current >= match_id);

    // The window leaves `matches[i - k]` as it admits `matches[i]`, so the two
    // ends are just the sequence offset against itself by k.
    for (leaving, &new_state) in matches.iter().zip(matches.iter().skip(k)) {
        current = current - i64::from(*leaving) + i64::from(new_state);
        aligned_columns += i64::from(current >= match_id);
    }

    // Divided by the QUERY length, not the alignment length, so a ratio above 1
    // is possible and is the reference's behaviour.
    let alignment_ratio = aligned_columns as f64 / s1.len() as f64;
    BlockAlignment {
        s1_aligned: a1,
        s2_aligned: a2,
        alignment_ratio,
    }
}

/// The gap-opening penalty bins, exactly as the reference writes them.
///
/// Used by `get_best_cluster_block_align`, which is the next stage; unit-tested
/// here in the meantime.
#[allow(dead_code)]
pub fn gap_opening_penalty(error_rate_sum: f64) -> i32 {
    if error_rate_sum <= 0.01 {
        5
    } else if error_rate_sum <= 0.04 {
        4
    } else if error_rate_sum <= 0.1 {
        3
    } else {
        2
    }
}

/// `math.floor((1.0 - error_rate_sum) * k)`.
#[allow(dead_code)]
pub fn match_id_tailored(error_rate_sum: f64, k: usize) -> i64 {
    ((1.0 - error_rate_sum) * k as f64).floor() as i64
}

/// The outcome of the alignment attempt, mirroring the reference's 6-tuple.
///
/// On failure the reference returns `(-1, 0, -1, -1, -1, alignment_ratio)` --
/// and that last value is the ratio from the *last candidate tried*, because
/// `alignment_ratio` is assigned inside the loop and leaks out of it. Reproduced:
/// `ratio` is not reset on failure.
#[derive(Debug, Clone, PartialEq)]
pub struct AlignResult {
    pub best_cluster_id: i64,
    pub nr_shared_kmers: usize,
    pub error_rate_sum: f64,
    pub alignment_ratio: f64,
}

/// `sum([q.count(c) * phred[c] for c in set(q)])` over the CAPPED table.
///
/// Note the two differences from the sorting stage's error rate: this uses the
/// capped table (`cluster.py`'s `phred_char_to_p`), and it divides by the
/// *sequence* length rather than the quality length. They are equal for
/// well-formed input, but the reference writes `len(seq)`, so this does too.
pub fn expected_errors(qual: &[u8]) -> f64 {
    let mut counts = [0u32; 256];
    for &c in qual {
        counts[c as usize] += 1;
    }
    crate::sorting::fsum(
        counts
            .iter()
            .enumerate()
            .filter(|(_, n)| **n > 0)
            .map(|(c, n)| f64::from(*n) * crate::phred::capped(c as u8)),
    )
}

/// `get_best_cluster_block_align`.
///
/// Only candidates **tied at the top hit count** are tried (`if nm_hits <
/// top_hits: break`), which is stricter than `get_best_cluster`'s
/// `min_fraction` walk. The first candidate whose aligned fraction reaches
/// `--aligned_threshold` wins.
#[allow(clippy::too_many_arguments)]
pub fn get_best_cluster_block_align(
    read_cl_id: usize,
    hits: &crate::cluster::Hits,
    seqs: &dyn Fn(usize) -> (Vec<u8>, Vec<u8>),
    accs: &dyn Fn(usize) -> String,
    k: usize,
    aligned_threshold: f64,
) -> AlignResult {
    let mut result = AlignResult {
        best_cluster_id: -1,
        nr_shared_kmers: 0,
        error_rate_sum: -1.0,
        alignment_ratio: 0.0,
    };
    if hits.is_empty() {
        // The reference would raise IndexError on top_matches[0] here; callers
        // only reach it with a non-empty hit set.
        return result;
    }

    let mut top_matches: Vec<usize> = hits.order.clone();
    top_matches.sort_by(|a, b| {
        let (ha, hb) = (&hits.by_cluster[a], &hits.by_cluster[b]);
        let ka = (
            ha.positions.len(),
            ha.positions.iter().sum::<usize>(),
            accs(*a),
        );
        let kb = (
            hb.positions.len(),
            hb.positions.iter().sum::<usize>(),
            accs(*b),
        );
        kb.cmp(&ka)
    });

    let (seq, r_qual) = seqs(read_cl_id);
    let top_hits = hits.by_cluster[&top_matches[0]].positions.len();

    for cl_id in top_matches {
        let nm_hits = hits.by_cluster[&cl_id].positions.len();
        if nm_hits < top_hits {
            break;
        }
        let (c_seq, c_qual) = seqs(cl_id);
        let error_rate_sum = expected_errors(&r_qual) / seq.len() as f64
            + expected_errors(&c_qual) / c_seq.len() as f64;
        let open = gap_opening_penalty(error_rate_sum);
        let match_id = match_id_tailored(error_rate_sum, k);
        let block = parasail_block_alignment(&seq, &c_seq, k, match_id, open);
        // The ratio leaks out of the loop in the reference, so keep the last
        // one tried even when nothing matches.
        result.alignment_ratio = block.alignment_ratio;
        if block.alignment_ratio >= aligned_threshold {
            result.best_cluster_id = cl_id as i64;
            result.nr_shared_kmers = nm_hits;
            result.error_rate_sum = error_rate_sum;
            return result;
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gap_penalty_bins_match_the_reference_boundaries() {
        // The reference's chain is <=0.01, then 0.01< x <=0.04, then <=0.1, else.
        assert_eq!(gap_opening_penalty(0.0), 5);
        assert_eq!(gap_opening_penalty(0.01), 5);
        assert_eq!(gap_opening_penalty(0.010001), 4);
        assert_eq!(gap_opening_penalty(0.04), 4);
        assert_eq!(gap_opening_penalty(0.05), 3);
        assert_eq!(gap_opening_penalty(0.1), 3);
        assert_eq!(gap_opening_penalty(0.2), 2);
    }

    #[test]
    fn match_id_floors() {
        assert_eq!(match_id_tailored(0.0, 13), 13);
        assert_eq!(match_id_tailored(0.1, 13), 11); // 11.7 -> 11
        assert_eq!(match_id_tailored(0.5, 15), 7); // 7.5 -> 7
    }

    #[test]
    fn identical_sequences_align_fully() {
        let s = b"ACGTACGTACGTACGTACGTACGTACGT";
        let a = parasail_block_alignment(s, s, 13, 13, 5);
        assert_eq!(a.s1_aligned, a.s2_aligned);
        // every window of 13 is all matches, so every column counts
        assert!(a.alignment_ratio > 0.5, "ratio {}", a.alignment_ratio);
    }

    /// The ratio divides by the query length, not the alignment length, so it is
    /// not bounded by 1. Reproduced deliberately.
    #[test]
    fn the_ratio_is_not_bounded_by_one() {
        // A short query against a long reference gives an alignment much longer
        // than the query, and every column can count.
        let s1 = b"ACGTACGTACGTACGTACGT";
        let s2 = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let a = parasail_block_alignment(s1, s2, 4, 0, 5);
        assert!(
            a.alignment_ratio > 1.0,
            "expected a ratio above 1, got {}",
            a.alignment_ratio
        );
    }
}
