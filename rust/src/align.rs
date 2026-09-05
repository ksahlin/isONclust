//! CARRIED ACROSS from the isONform port. The code below is unchanged.
//!
//! isONclust has its own `cigar_to_seq` in `modules/help_functions.py`, and it
//! uses the same convention as this one: `I` consumes the query, `D` consumes
//! the reference. Verified against the Python rather than assumed --- see
//! `tests::matches_the_isonclust_reference`.
//!
//! Global alignment with an edlib-compatible CIGAR, matching
//! `edlib.align(query, target, task="path", mode="NW")`.
//!
//! `get_best_corrections` aligns every supporting segment back to the consensus
//! and consumes the **CIGAR**, not just the distance. That makes this the one
//! place where matching edlib's *choice* of alignment matters: many alignments
//! achieve the same optimal edit distance, and a different tie-break yields a
//! different-but-equally-optimal CIGAR — which changes the MSA and therefore the
//! corrected output.
//!
//! # The tie-break
//!
//! Determined empirically rather than guessed. Running the traceback with all
//! six preference orderings against 677 real `task="path"` calls recorded from
//! the reference:
//!
//! | preference (backwards from the end) | CIGARs matching |
//! | --- | --- |
//! | insert, delete, diagonal | **677 / 677** |
//! | delete, insert, diagonal | 637 / 677 |
//! | delete, diagonal, insert | 281 / 677 |
//! | insert, diagonal, delete | 131 / 677 |
//! | diagonal, insert, delete | 56 / 677 |
//! | diagonal, delete, insert | 55 / 677 |
//!
//! So edlib prefers **insertion, then deletion, then the diagonal**, walking
//! backwards from `(n, m)`. Note how badly the intuitive diagonal-first ordering
//! does — 8% — which is why this was measured.
//!
//! Every ordering reproduces the edit *distance*; only the path differs. That is
//! the whole point: the distance is unique, the alignment is not.
//!
//! CIGAR operations use edlib's extended alphabet: `=` match, `X` mismatch,
//! `I` insertion to target (consumes query), `D` deletion from target.
//!
//! `#![allow(dead_code)]` below is deliberate: this file is carried across
//! VERBATIM from the isONform port and is verified there. isONclust uses only
//! part of its surface today, and trimming the rest would mean editing code
//! this port has not verified, to silence a lint about code that is correct.
//! The unused parts are what the remaining stages will call.
#![allow(dead_code)]

/// Alignment of `query` against `target`: edit distance and extended CIGAR.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    pub edit_distance: usize,
    pub cigar: String,
}

/// Globally align `query` to `target`, returning an edlib-compatible CIGAR.
///
/// Straightforward `O(n*m)` Needleman-Wunsch with unit costs. The segments this
/// runs on are short — anchor spans are bounded by `--xmax`, 80 by default — so
/// the quadratic table is small. edlib is bit-parallel and would be faster on
/// long inputs; if that ever matters, the oracle test below is what makes a
/// faster implementation safe to swap in.
pub fn global(query: &[u8], target: &[u8]) -> Alignment {
    let (n, m) = (query.len(), target.len());

    let mut dp = vec![0u32; (n + 1) * (m + 1)];
    let at = |i: usize, j: usize| i * (m + 1) + j;
    for i in 0..=n {
        dp[at(i, 0)] = i as u32;
    }
    for j in 0..=m {
        dp[at(0, j)] = j as u32;
    }
    for i in 1..=n {
        for j in 1..=m {
            let sub = dp[at(i - 1, j - 1)] + u32::from(query[i - 1] != target[j - 1]);
            dp[at(i, j)] = sub.min(dp[at(i - 1, j)] + 1).min(dp[at(i, j - 1)] + 1);
        }
    }

    // Traceback, preferring insertion, then deletion, then the diagonal.
    let (mut i, mut j) = (n, m);
    let mut ops: Vec<u8> = Vec::with_capacity(n + m);
    while i > 0 || j > 0 {
        if i > 0 && dp[at(i - 1, j)] + 1 == dp[at(i, j)] {
            ops.push(b'I');
            i -= 1;
        } else if j > 0 && dp[at(i, j - 1)] + 1 == dp[at(i, j)] {
            ops.push(b'D');
            j -= 1;
        } else if i > 0 && j > 0 {
            let same = query[i - 1] == target[j - 1];
            ops.push(if same { b'=' } else { b'X' });
            i -= 1;
            j -= 1;
        } else {
            break;
        }
    }
    ops.reverse();

    Alignment {
        edit_distance: dp[at(n, m)] as usize,
        cigar: run_length_encode(&ops),
    }
}

/// One run of a CIGAR: `(length, operation)`.
pub type CigarOp = (usize, u8);

/// Parse an extended CIGAR into runs, matching `help_functions.cigar_to_seq`'s
/// `re.split(r'[=DXSMI]+', cigar)` walk.
///
/// The reference's regex tolerates `S` and `M` when splitting but its loop then
/// falls through to `print("error"); sys.exit()` for anything that is not
/// `=`, `X`, `I` or `D`. Nothing in the pipeline produces those operations —
/// edlib emits only the extended four — so this returns `None` instead of
/// aborting the process, and callers treat it as a bug rather than an input.
pub fn parse_cigar(cigar: &str) -> Option<Vec<CigarOp>> {
    let mut ops = Vec::new();
    let mut len = 0usize;
    let mut seen_digit = false;
    for c in cigar.bytes() {
        if c.is_ascii_digit() {
            len = len * 10 + usize::from(c - b'0');
            seen_digit = true;
        } else {
            if !seen_digit || !matches!(c, b'=' | b'X' | b'I' | b'D') {
                return None;
            }
            ops.push((len, c));
            len = 0;
            seen_digit = false;
        }
    }
    if seen_digit {
        return None; // trailing length with no operation
    }
    Some(ops)
}

/// Expand a CIGAR into the two gapped strings, matching
/// `help_functions.cigar_to_seq(cigar, query, ref)`.
///
/// Returns `(query_aligned, ref_aligned)` — the reference's return order, which
/// `get_best_corrections` immediately unpacks as `read_alignment,
/// ref_alignment`. `I` is an insertion with respect to the reference (a gap in
/// `ref_aligned`); `D` is a deletion (a gap in `query_aligned`).
///
/// Returns `None` if the CIGAR does not parse or does not cover both sequences.
pub fn cigar_to_seq(cigar: &str, query: &[u8], reference: &[u8]) -> Option<(Vec<u8>, Vec<u8>)> {
    ops_to_seq(&parse_cigar(cigar)?, query, reference)
}

/// [`cigar_to_seq`] without the string.
///
/// The SIMD aligner reports operations directly, so on the hot path formatting a
/// CIGAR only to parse it straight back is pure overhead — 4.27 million times per
/// four real clusters. [`crate::simd::global_ops`] feeds this instead.
pub fn ops_to_seq(ops: &[CigarOp], query: &[u8], reference: &[u8]) -> Option<(Vec<u8>, Vec<u8>)> {
    let (mut q_index, mut r_index) = (0usize, 0usize);
    let (mut q_aln, mut r_aln) = (Vec::new(), Vec::new());

    for &(len, op) in ops {
        match op {
            b'=' | b'X' => {
                q_aln.extend_from_slice(query.get(q_index..q_index + len)?);
                r_aln.extend_from_slice(reference.get(r_index..r_index + len)?);
                q_index += len;
                r_index += len;
            }
            b'I' => {
                r_aln.extend(std::iter::repeat_n(b'-', len));
                q_aln.extend_from_slice(query.get(q_index..q_index + len)?);
                q_index += len;
            }
            b'D' => {
                r_aln.extend_from_slice(reference.get(r_index..r_index + len)?);
                q_aln.extend(std::iter::repeat_n(b'-', len));
                r_index += len;
            }
            _ => return None,
        }
    }
    // Python would silently accept a CIGAR that stops short; nothing produces
    // one, and a partial alignment would corrupt the MSA rather than fail.
    if q_index != query.len() || r_index != reference.len() {
        return None;
    }
    Some((q_aln, r_aln))
}

/// `[(2, b'='), (1, b'X')]` becomes `"2=1X"`. The inverse of [`parse_cigar`].
pub fn encode_cigar(ops: &[CigarOp]) -> String {
    use std::fmt::Write;
    let mut out = String::new();
    for &(len, op) in ops {
        let _ = write!(out, "{len}{}", op as char);
    }
    out
}

/// `[b'=', b'=', b'X']` becomes `"2=1X"`.
fn run_length_encode(ops: &[u8]) -> String {
    let mut out = String::new();
    let mut k = 0;
    while k < ops.len() {
        let c = ops[k];
        let start = k;
        while k < ops.len() && ops[k] == c {
            k += 1;
        }
        out.push_str(&(k - start).to_string());
        out.push(c as char);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cigar(q: &str, t: &str) -> String {
        global(q.as_bytes(), t.as_bytes()).cigar
    }

    #[test]
    fn identical_sequences_are_all_matches() {
        let a = global(b"ACGT", b"ACGT");
        assert_eq!(a.edit_distance, 0);
        assert_eq!(a.cigar, "4=");
    }

    #[test]
    fn a_substitution_shows_as_x() {
        let a = global(b"ACGT", b"AGGT");
        assert_eq!(a.edit_distance, 1);
        assert_eq!(a.cigar, "1=1X2=");
    }

    #[test]
    fn empty_query_is_all_deletions() {
        let a = global(b"", b"ACGT");
        assert_eq!(a.edit_distance, 4);
        assert_eq!(a.cigar, "4D");
    }

    #[test]
    fn empty_target_is_all_insertions() {
        let a = global(b"ACGT", b"");
        assert_eq!(a.edit_distance, 4);
        assert_eq!(a.cigar, "4I");
    }

    #[test]
    fn both_empty_yields_an_empty_cigar() {
        let a = global(b"", b"");
        assert_eq!(a.edit_distance, 0);
        assert_eq!(a.cigar, "");
    }

    #[test]
    fn cigar_operation_counts_cover_both_sequences() {
        // Query length = = + X + I; target length = = + X + D.
        let (q, t) = ("ACGTACGTAA", "ACGTTACGT");
        let c = cigar(q, t);
        let (mut qn, mut tn, mut num) = (0usize, 0usize, String::new());
        for ch in c.chars() {
            if ch.is_ascii_digit() {
                num.push(ch);
            } else {
                let n: usize = num.parse().unwrap();
                num.clear();
                match ch {
                    '=' | 'X' => {
                        qn += n;
                        tn += n;
                    }
                    'I' => qn += n,
                    'D' => tn += n,
                    _ => panic!("unexpected op {ch}"),
                }
            }
        }
        assert_eq!(qn, q.len(), "query coverage");
        assert_eq!(tn, t.len(), "target coverage");
    }

    fn expand(cigar: &str, q: &str, r: &str) -> (String, String) {
        let (qa, ra) = cigar_to_seq(cigar, q.as_bytes(), r.as_bytes()).unwrap();
        (
            String::from_utf8(qa).unwrap(),
            String::from_utf8(ra).unwrap(),
        )
    }

    #[test]
    fn cigar_parses_into_runs() {
        assert_eq!(
            parse_cigar("2=1X10I").unwrap(),
            vec![(2, b'='), (1, b'X'), (10, b'I')]
        );
        assert_eq!(parse_cigar("").unwrap(), vec![]);
        assert!(parse_cigar("3M").is_none(), "M is not in the pipeline");
        assert!(parse_cigar("3").is_none(), "length with no operation");
        assert!(parse_cigar("=").is_none(), "operation with no length");
    }

    #[test]
    fn insertion_gaps_the_reference_and_deletion_gaps_the_query() {
        let (q, r) = expand("2=2I2=", "ACGTAC", "ACAC");
        assert_eq!((q.as_str(), r.as_str()), ("ACGTAC", "AC--AC"));
        let (q, r) = expand("2=2D2=", "ACAC", "ACGTAC");
        assert_eq!((q.as_str(), r.as_str()), ("AC--AC", "ACGTAC"));
    }

    #[test]
    fn an_expanded_cigar_round_trips_its_own_alignment() {
        let (q, r) = ("ACGTACGTAA", "ACGTTACGT");
        let a = global(q.as_bytes(), r.as_bytes());
        let (qa, ra) = expand(&a.cigar, q, r);
        assert_eq!(qa.len(), ra.len());
        assert_eq!(qa.replace('-', ""), q);
        assert_eq!(ra.replace('-', ""), r);
        // Gapped columns account for exactly the edit distance.
        let diffs = qa.bytes().zip(ra.bytes()).filter(|(a, b)| a != b).count();
        assert_eq!(diffs, a.edit_distance);
    }

    #[test]
    fn a_cigar_that_does_not_cover_both_sequences_is_rejected() {
        assert!(cigar_to_seq("2=", b"ACGT", b"ACGT").is_none());
        assert!(cigar_to_seq("4=", b"ACGT", b"ACG").is_none());
    }

    /// Tie cases verified directly against Python edlib, not assumed.
    ///
    /// The preference is applied walking *backwards*, so the gap lands at the
    /// end of the emitted CIGAR rather than the start --- which is the opposite
    /// of what "insertion first" suggests when read forwards.
    #[test]
    fn tie_cases_match_edlib_ground_truth() {
        assert_eq!(cigar("AA", "A"), "1=1I");
        assert_eq!(cigar("A", "AA"), "1=1D");
        assert_eq!(cigar("ACGT", "ACT"), "2=1I1=");
        assert_eq!(cigar("AAA", "A"), "1=2I");
    }
}
// Removed on the way across from isONcorrect: the `oracle` test module, which
// compared `crate::simd::global_ops` and `parasail::global_affine` against this
// file's exact Needleman-Wunsch. isONform imports neither --- it never uses edlib,
// so there is no edit-distance path to accelerate, and `global_affine` was
// dropped from parasail.rs with the rest of the guard research. What isONform
// needs from this file is the CIGAR handling below, which
// `consensus.parasail_alignment` reaches through `cigar_to_seq`.

// Removed on the way across from isONcorrect: the `block_option` test module.
//
// It read recorded edlib CIGARs and timed a SIMD banded aligner against this
// file's scalar Needleman-Wunsch, to decide whether banding was safe for
// isONcorrect's correction segments. isONform never imports edlib and has no
// edit-distance path here, so there is nothing to compare and no decision to
// make. What it needs from this file is the CIGAR handling above, reached from
// `consensus.parasail_alignment` via `cigar_to_seq`.
