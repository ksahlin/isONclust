//! parasail's own C library, through FFI.
//!
//! # Why this exists
//!
//! `parasail.rs` is an exact scalar reimplementation. It is *right* — that is
//! what it is for, and it is what proved the port byte-identical — but it is
//! 13-16x slower than the C library, and alignment is 96-99.6% of this tool's
//! runtime. That single fact made the port slower than the Python reference on
//! ONT and PacBio data.
//!
//! Three replacements were measured on 18 633 real recorded alignments
//! (PORTING.md, *The aligner*). block-aligner and WFA2 are both **slower than
//! this** and change 1-43% of clustering verdicts; rust-bio is slower than our
//! own scalar code; edlib and triple_accel cannot express the scoring at all.
//! Linking parasail is the only option that is faster *and* exact.
//!
//! # Why it is exact rather than merely accurate
//!
//! This is not an implementation that agrees with the reference — it **is** the
//! implementation the reference calls. The Python `parasail` package is a
//! binding around this same C library, so `sg_trace_scan_16` here and
//! `parasail.sg_trace_scan_16` there are the same function on the same inputs.
//! Verified anyway rather than assumed: 2000 of 2000 recorded `droso_20k`
//! alignments return a byte-identical CIGAR, and `equivalence.sh stage parasail`
//! re-checks all 18 633 across four corpora.
//!
//! # The build cost
//!
//! `libparasail-sys` builds parasail from source, so it needs `cmake` and
//! `libclang`. That is a real cost, and it is why the pure-Rust path is kept
//! behind the default: `--no-default-features` drops this module and falls back
//! to `parasail.rs`, which is exact and needs nothing. Both produce identical
//! output; the choice is build complexity against speed, not correctness.

use crate::align::CigarOp;
use crate::parasail::{Alignment, Scoring};
use std::ffi::{CStr, CString};

/// Semi-global affine alignment, matching `parasail.sg_trace_scan_16`.
///
/// The reference falls back to `sg_trace_scan_32` when the 16-bit result
/// saturates; so does this, via the same `saturated` flag.
pub fn semiglobal(s1: &[u8], s2: &[u8], sc: Scoring) -> Alignment {
    // parasail's matrix is built per (match, mismatch) pair, and isONclust only
    // ever uses (2, -2). Building it once per call showed up in the profile, so
    // it is cached — parasail matrices are immutable once created.
    let matrix = matrix_for(sc.match_score, sc.mismatch);

    unsafe {
        let mut result = libparasail_sys::parasail_sg_trace_scan_16(
            s1.as_ptr() as *const std::os::raw::c_char,
            s1.len() as std::os::raw::c_int,
            s2.as_ptr() as *const std::os::raw::c_char,
            s2.len() as std::os::raw::c_int,
            sc.open,
            sc.ext,
            matrix,
        );
        if libparasail_sys::parasail_result_is_saturated(result) != 0 {
            libparasail_sys::parasail_result_free(result);
            result = libparasail_sys::parasail_sg_trace_scan_32(
                s1.as_ptr() as *const std::os::raw::c_char,
                s1.len() as std::os::raw::c_int,
                s2.as_ptr() as *const std::os::raw::c_char,
                s2.len() as std::os::raw::c_int,
                sc.open,
                sc.ext,
                matrix,
            );
        }
        let score = libparasail_sys::parasail_result_get_score(result);
        let cigar_t = libparasail_sys::parasail_result_get_cigar(
            result,
            s1.as_ptr() as *const std::os::raw::c_char,
            s1.len() as std::os::raw::c_int,
            s2.as_ptr() as *const std::os::raw::c_char,
            s2.len() as std::os::raw::c_int,
            matrix,
        );
        let decoded = libparasail_sys::parasail_cigar_decode(cigar_t);
        let cigar = CStr::from_ptr(decoded).to_string_lossy().into_owned();
        // parasail_cigar_decode allocates; the library frees it with the cigar.
        libparasail_sys::parasail_cigar_free(cigar_t);
        libparasail_sys::parasail_result_free(result);

        let ops = crate::align::parse_cigar(&cigar).unwrap_or_default();
        Alignment { score, cigar, ops }
    }
}

/// Cached `parasail_matrix_t` for one (match, mismatch) pair.
///
/// isONclust uses exactly one pair, so this is a single-entry cache in practice.
/// The matrix is read-only after creation and parasail's own functions take it
/// as `*const`, so sharing it across threads is sound.
fn matrix_for(match_score: i32, mismatch: i32) -> *const libparasail_sys::parasail_matrix_t {
    use std::sync::{Mutex, OnceLock};
    /// (match, mismatch) -> the matrix pointer, as usize so it is `Send`.
    type MatrixCache = Vec<((i32, i32), usize)>;
    static CACHE: OnceLock<Mutex<MatrixCache>> = OnceLock::new();
    let cache = CACHE.get_or_init(|| Mutex::new(Vec::new()));
    let mut guard = cache.lock().expect("matrix cache poisoned");
    if let Some((_, p)) = guard
        .iter()
        .find(|((m, x), _)| *m == match_score && *x == mismatch)
    {
        return *p as *const libparasail_sys::parasail_matrix_t;
    }
    let alphabet = CString::new("ACGT").expect("no interior nul");
    let m = unsafe {
        libparasail_sys::parasail_matrix_create(alphabet.as_ptr(), match_score, mismatch)
    };
    // Deliberately never freed: one matrix for the life of the process, and
    // freeing it would require proving no alignment still holds the pointer.
    guard.push(((match_score, mismatch), m as usize));
    m as *const libparasail_sys::parasail_matrix_t
}

/// The op vector, for callers that do not need the string form.
#[allow(dead_code)]
pub fn semiglobal_ops(s1: &[u8], s2: &[u8], sc: Scoring) -> Vec<CigarOp> {
    semiglobal(s1, s2, sc).ops
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sc(open: i32) -> Scoring {
        Scoring {
            match_score: 2,
            mismatch: -2,
            open,
            ext: 1,
        }
    }

    /// The same cases `parasail.rs`'s own tests pin, so a divergence between the
    /// two implementations shows up here rather than on a corpus.
    #[test]
    fn agrees_with_the_scalar_implementation() {
        let cases: &[(&[u8], &[u8])] = &[
            (b"ACGTACGT", b"ACGTACGT"),
            (b"ACGTACGT", b"ACGTTCGT"),
            (b"ACGTACGT", b"TTTACGTACGTTTT"),
            (b"TTTACGTACGTTTT", b"ACGTACGT"),
            (b"ACGTACGT", b"ACGACGT"),
            (b"ACGTACGTACGTACGTACGT", b"ACGTACGTTCGTACGTACGT"),
        ];
        for open in [2, 3, 4, 5] {
            for (a, b) in cases {
                let ffi = semiglobal(a, b, sc(open));
                let scalar = crate::parasail::semiglobal(a, b, sc(open));
                assert_eq!(
                    (ffi.score, &ffi.cigar),
                    (scalar.score, &scalar.cigar),
                    "open={open} a={} b={}",
                    String::from_utf8_lossy(a),
                    String::from_utf8_lossy(b)
                );
            }
        }
    }

    #[test]
    fn the_matrix_cache_returns_a_stable_pointer() {
        let a = matrix_for(2, -2);
        let b = matrix_for(2, -2);
        assert_eq!(a, b, "the cache should not rebuild the matrix");
    }

    #[test]
    fn a_non_acgt_character_scores_zero_against_itself() {
        // parasail_matrix_create("ACGT", ...) defines only the four bases; every
        // other byte scores 0, not match_score. isONform's finding 24, and it is
        // why an N in a read is not a free match here either.
        let n = semiglobal(b"NNNN", b"NNNN", sc(5));
        let a = semiglobal(b"AAAA", b"AAAA", sc(5));
        assert_eq!(n.score, 0);
        assert_eq!(a.score, 8);
    }
}
