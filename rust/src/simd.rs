//! The SIMD aligner: block-aligner behind a parasail-shaped call.
//!
//! # Why
//!
//! Profiling says alignment is the dominant cost of the whole tool, in *two*
//! stages, and the graph is not the problem:
//!
//! | stage | 1 000-read batch | 250-read batch, 7 343-node graph |
//! | --- | --- | --- |
//! | intervals | 2.33s | 0.25s |
//! | graph build | **0.35s** | **0.03s** |
//! | simplify (`pop_bubbles`) | 2.29s | **46.04s** |
//! | merge (`merge_consensuses`) | **58.05s** | 20.56s |
//!
//! Which of simplify and merge dominates depends on the graph, but sampling both
//! lands in `parasail::semiglobal_with` and `parasail::traceback`. So one aligner
//! swap addresses the whole pipeline.
//!
//! isONcorrect already made this swap and measured **8-99x**, on 1 400 verified
//! alignments. This port has 54 884 recorded real-parasail cases to check against,
//! which is 39x the evidence.
//!
//! # Free end gaps: what parasail does, and what block-aligner gives
//!
//! [`crate::parasail::semiglobal`] zeroes row 0 **and** column 0, and its end-cell
//! scan looks at the last row **and** the last column. Both sequences may dangle
//! at both ends --- an overlap alignment. That is what makes `align_to_merge`
//! work on isoforms whose lengths differ.
//!
//! block-aligner 0.5.1 parameterises `Block` as
//! `Block<TRACE, X_DROP, LOCAL_START, FREE_QUERY_START_GAPS, FREE_QUERY_END_GAPS>`
//! --- free gaps on the **query** only, not symmetric.
//!
//! With the query as the *shorter* sequence, freeing the query's start and end
//! gaps lets the shorter sequence sit anywhere inside the longer one, which is the
//! case `align_to_merge` cares about: it always compares a long consensus against
//! a short one and asks whether the short is contained. The asymmetric form is
//! therefore the right shape, *provided* the caller keeps the orientation. This
//! module enforces it rather than trusting callers --- see [`semiglobal_ops`].
//!
//! What it is **not** is bit-identical to parasail. Two reasons, and both are
//! measured rather than assumed by [`oracle`]:
//!
//! * the band is adaptive, so a score can in principle be sub-optimal;
//! * even at an equal score, a different equally-optimal path may be reported,
//!   which changes the CIGAR and therefore the diversity computation that decides
//!   a merge.
//!
//! Reference identity is no longer the goal, so a different-but-equally-optimal
//! path is acceptable. A *worse* score is not, and that is what the oracle gates.
//!
//! CARRIED ACROSS from the isONform port, VERBATIM apart from this header and
//! the removal of a submodule that referenced isONform's own call sites.
//! `#![allow(dead_code)]` for the same reason as `parasail.rs`: isONclust
//! exercises only part of the surface, and trimming the rest would mean editing
//! verified code to silence a lint about code that is correct.
//!
//! This module is EXPERIMENTAL here and off by default. `PORTING.md`'s aligner
//! section has the measurements: on isONclust's own recorded alignments it is
//! slower than parasail's C library and changes 9-43% of clustering verdicts, so
//! it is kept for the record rather than as a candidate.
#![allow(dead_code)]

use crate::align::CigarOp;
use crate::parasail::Scoring;

/// Upper bound on the block size, and therefore on the query length this can
/// handle. Block sizes must be powers of two and below 2^16.
///
/// # Free end gaps forbid banding
///
/// block-aligner asserts `min_size > query.len()` whenever
/// `FREE_QUERY_END_GAPS` is set, so the block must span the whole query and there
/// is no adaptive band. That is a real cost: isONcorrect's 8-99x came from SIMD
/// *and* banding on short segments, and only the SIMD half is available here.
/// Measured rather than assumed --- see `oracle` and the benchmark.
const MAX_BLOCK: usize = 32768;
/// Starting capacity and growth granularity, so a run of slightly longer
/// sequences does not reallocate every call.
const CAP_STEP: usize = 512;

/// The block size for a given query length: the smallest power of two strictly
/// greater than it, which is what the `FREE_QUERY_END_GAPS` assert demands.
fn block_size_for(query_len: usize) -> usize {
    (query_len + 1).next_power_of_two().clamp(32, MAX_BLOCK)
}

/// Semi-global alignment ops for a pair of sequences, SIMD.
///
/// **Orientation is enforced, not requested.** The shorter sequence becomes the
/// query, because block-aligner's free-gap flags apply to the query only. When the
/// arguments arrive the other way round the ops are transposed on the way out, so
/// callers see ops describing `(s1, s2)` in the order they passed them --- an
/// insertion stays an insertion relative to *their* first argument.
///
/// Returns `false` when the pair cannot be aligned this way (empty input, a score
/// that will not fit `i8`, or a truncated result), and the caller should fall back
/// to the exact path. Silently returning a wrong alignment is the one thing that
/// must not happen, since a merge decision hangs on it.
pub fn semiglobal_ops(s1: &[u8], s2: &[u8], sc: Scoring, out: &mut Vec<CigarOp>) -> Option<i32> {
    out.clear();
    if s1.is_empty() || s2.is_empty() {
        return None;
    }
    let swapped = s1.len() < s2.len();
    let (query, reference) = if swapped { (s1, s2) } else { (s2, s1) };
    let score = SCRATCH.with_borrow_mut(|s| s.align_into(query, reference, sc, out))?;
    // `out` describes (query, reference) = (shorter, longer). The caller passed
    // (s1, s2). Transpose when those disagree: an insertion in one orientation is
    // a deletion in the other.
    let caller_query_is_shorter = s1.len() <= s2.len();
    let produced_query_is_s1 = swapped;
    if produced_query_is_s1 != caller_query_is_shorter || !produced_query_is_s1 {
        for op in out.iter_mut() {
            op.1 = match op.1 {
                b'I' => b'D',
                b'D' => b'I',
                other => other,
            };
        }
    }
    Some(score)
}

/// `NucMatrix` with parasail's catch-all row: any non-ACGT byte scores 0.
///
/// `NucMatrix`'s alphabet is `A T C G N`, and `set` writes both `(a,b)` and
/// `(b,a)`, so zeroing the `N` row and column is enough. Bytes outside the
/// alphabet cannot appear --- `convert_char` asserts `A..=Z`, and the caller
/// refuses anything else.
fn nuc_matrix_with_parasail_catchall(
    match_score: i8,
    mismatch: i8,
) -> block_aligner::scores::NucMatrix {
    use block_aligner::scores::{Matrix, NucMatrix};
    let mut m = NucMatrix::new_simple(match_score, mismatch);
    for b in *b"ATCGN" {
        m.set(b'N', b, 0);
    }
    m
}

thread_local! {
    /// Reusable allocations. **Not a micro-optimisation:** the exact path
    /// allocates and zeroes an `n * m` byte table per call --- ~1.9 MB at the
    /// measured mean consensus length of 1 364 bp --- and there are thousands of
    /// calls per batch. isONcorrect measured per-call setup at 18.84s against
    /// 11.54s reused on the same work.
    static SCRATCH: std::cell::RefCell<Scratch> = std::cell::RefCell::new(Scratch::new());
}

struct Scratch {
    /// `<TRACE, X_DROP, LOCAL_START, FREE_QUERY_START_GAPS, FREE_QUERY_END_GAPS>`
    /// --- trace on (the CIGAR is needed), and the query's end gaps free, which is
    /// the semi-global shape parasail's `sg` provides. See the module docs.
    block: block_aligner::scan_block::Block<true, false, false, true, true>,
    q: block_aligner::scan_block::PaddedBytes,
    r: block_aligner::scan_block::PaddedBytes,
    cigar: block_aligner::cigar::Cigar,
    cap: usize,
    /// The block size the current allocations were built for.
    blk: usize,
}

impl Scratch {
    fn new() -> Self {
        use block_aligner::{scan_block::*, scores::*};
        let blk = block_size_for(CAP_STEP);
        Self {
            block: Block::new(CAP_STEP, CAP_STEP, blk),
            q: PaddedBytes::new::<NucMatrix>(CAP_STEP, blk),
            r: PaddedBytes::new::<NucMatrix>(CAP_STEP, blk),
            cigar: block_aligner::cigar::Cigar::new(CAP_STEP, CAP_STEP),
            cap: CAP_STEP,
            blk,
        }
    }

    /// `Block::new` takes *maximum* lengths and `align` the actual ones, so one
    /// block sized to the largest pair seen serves every smaller one. Grows only.
    fn reserve(&mut self, len: usize, query_len: usize) -> usize {
        use block_aligner::{scan_block::*, scores::*};
        let want_blk = block_size_for(query_len);
        if len <= self.cap && want_blk <= self.blk {
            return self.blk;
        }
        self.cap = len.max(self.cap).next_multiple_of(CAP_STEP);
        self.blk = want_blk.max(self.blk);
        self.block = Block::new(self.cap, self.cap, self.blk);
        self.q = PaddedBytes::new::<NucMatrix>(self.cap, self.blk);
        self.r = PaddedBytes::new::<NucMatrix>(self.cap, self.blk);
        self.cigar = block_aligner::cigar::Cigar::new(self.cap, self.cap);
        self.blk
    }

    fn align_into(
        &mut self,
        query: &[u8],
        reference: &[u8],
        sc: Scoring,
        out: &mut Vec<CigarOp>,
    ) -> Option<i32> {
        use block_aligner::{cigar::Operation, scores::*};

        // block-aligner takes i8 scores, and its gap penalties are negative.
        let (Ok(m), Ok(x), Ok(open), Ok(extend)) = (
            i8::try_from(sc.match_score),
            i8::try_from(sc.mismatch),
            i8::try_from(-sc.open),
            i8::try_from(-sc.ext),
        ) else {
            return None;
        };
        // parasail's `matrix_create("ACGT", m, x)` builds a 5x5 matrix whose
        // fifth row and column --- every byte outside ACGT --- is **zero**, not
        // the mismatch penalty (`PORTING.md` finding 24). `NucMatrix::new_simple`
        // scores `N` as a plain mismatch, and that single difference was the
        // whole of the earlier oracle failure: 20 000 cases scoring worse and
        // 8 812 higher, uniformly across every length bucket, because ambiguity
        // codes are scattered through the consensuses rather than concentrated by
        // length. Diagnosed by `diagnose_against_exact_parasail`, where the `N`
        // case was the only one of seven that disagreed.
        let matrix = nuc_matrix_with_parasail_catchall(m, x);
        let gaps = Gaps { open, extend };

        if query.len() >= MAX_BLOCK {
            return None; // beyond what a full-width block can cover
        }
        let blk = self.reserve(query.len().max(reference.len()), query.len());
        self.q.set_bytes::<NucMatrix>(query, blk);
        self.r.set_bytes::<NucMatrix>(reference, blk);
        // `min_size > query.len()` is required by FREE_QUERY_END_GAPS, so the
        // range is a single size: no band.
        let lo = block_size_for(query.len());
        self.block
            .align(&self.q, &self.r, &matrix, gaps, lo..=blk, 0);

        let res = self.block.res();
        // A result that did not consume the query is not the alignment asked for.
        if res.query_idx != query.len() {
            return None;
        }
        self.block.trace().cigar_eq(
            &self.q,
            &self.r,
            res.query_idx,
            res.reference_idx,
            &mut self.cigar,
        );

        out.reserve(self.cigar.len() + 4);
        for i in 0..self.cigar.len() {
            let ol = self.cigar.get(i);
            let op = match ol.op {
                Operation::Eq => b'=',
                Operation::X => b'X',
                Operation::I => b'I',
                Operation::D => b'D',
                Operation::M => b'M',
                _ => return None,
            };
            out.push((ol.len, op));
        }
        if out.is_empty() {
            return None;
        }

        // block-aligner's CIGAR covers only the aligned region: with free end
        // gaps the dangling prefix and suffix of the *reference* are simply not
        // in it. parasail's `sg` traceback does include them as I/D ops, and
        // `parse_cigar_diversity_isoform_level` needs them --- `mergeable_start`
        // and `mergeable_end` are computed from exactly those columns, which is
        // what `--delta_iso_len_5`/`_3` bound.
        //
        // So they are reconstructed here from the result indices rather than
        // left out. `res.query_idx`/`res.reference_idx` are the *end* positions,
        // and the CIGAR's own spans give the aligned extent, so the difference is
        // the leading dangle and the tail is whatever is left over.
        let (mut q_span, mut r_span) = (0usize, 0usize);
        for &(l, op) in out.iter() {
            match op {
                b'=' | b'X' | b'M' => {
                    q_span += l;
                    r_span += l;
                }
                b'I' => q_span += l,
                b'D' => r_span += l,
                _ => {}
            }
        }
        let q_lead = res.query_idx.saturating_sub(q_span);
        let r_lead = res.reference_idx.saturating_sub(r_span);
        let q_tail = query.len().saturating_sub(res.query_idx);
        let r_tail = reference.len().saturating_sub(res.reference_idx);
        if r_lead > 0 {
            out.insert(0, (r_lead, b'D'));
        }
        if q_lead > 0 {
            out.insert(0, (q_lead, b'I'));
        }
        if q_tail > 0 {
            out.push((q_tail, b'I'));
        }
        if r_tail > 0 {
            out.push((r_tail, b'D'));
        }
        Some(res.score)
    }
}

/// Does the SIMD aligner ever score *worse* than exact parasail?
///
/// The one question that matters. block-aligner's band is adaptive, so a
/// sub-optimal score is possible in principle; a different equally-optimal path
/// is expected and acceptable, since reference identity is no longer the goal.
/// So this gates the score and merely *reports* CIGAR differences.
///
/// ```text
/// PARASAIL_CASES=<dump>/parasail_cases.tsv \
///   cargo test --release --lib simd::oracle -- --nocapture
/// ```
///
/// isONcorrect validated the same swap on 1 400 alignments. This runs on 54 884.
#[cfg(test)]
pub mod oracle {
    use super::*;

    /// Does the SIMD aligner ever score *worse* than exact parasail?
    ///
    /// The one question that matters. block-aligner's band is adaptive, so a
    /// sub-optimal score is possible in principle; a different equally-optimal
    /// path is expected and acceptable now that reference identity is not the
    /// goal. So this gates the score and merely *reports* CIGAR differences.
    ///
    /// The recorded cases carry parasail's own score and scoring parameters, so
    /// the comparison is against what the real library returned rather than
    /// against a rescore of my own devising.
    ///
    /// isONcorrect validated the same swap on 1 400 alignments; this runs on
    /// 54 884.
    #[test]
    fn never_scores_worse_than_exact_parasail() {
        let Ok(path) = std::env::var("PARASAIL_CASES") else {
            eprintln!(
                "SKIPPED: set PARASAIL_CASES to a parasail_cases.tsv to check the \
                 SIMD aligner against the exact one. Without it this asserts nothing."
            );
            return;
        };
        let text = std::fs::read_to_string(&path).expect("cases readable");
        let (mut n, mut worse, mut better, mut refused) = (0u64, 0u64, 0u64, 0u64);
        let mut worst = 0i32;
        // Bucket by **identity**, not length. The merge only ever cares about
        // pairs similar enough to plausibly merge --- `align_to_merge` needs a
        // >=100bp shared region within the diversity budget --- so losing the
        // optimum on a dissimilar pair costs nothing. Identity comes from the
        // recorded CIGAR: matched columns over alignment columns.
        let mut buckets = [[0u64; 3]; 5]; // [identity bucket][worse, equal, better]
        let mut worse_by_delta = [0u64; 4]; // |delta| 1-2, 3-10, 11-50, 51+
        let mut ops = Vec::new();
        // Skip `#` comments and blanks, and **panic** on a malformed row rather
        // than defaulting. The first version of this used `unwrap_or` fallbacks,
        // so the header line parsed as data and every scoring parameter silently
        // became a plausible-looking default --- the mismatch is -8 in these
        // recordings, not the -2 the default assumed. The comparison then
        // measured nothing, and reported the same 20 000 failures before and
        // after a real fix. Same failure mode `ISONFORM_BUG_COMPAT` guards
        // against by exiting on an unknown name.
        for line in text
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let f: Vec<&str> = line.split('\t').collect();
            assert!(f.len() >= 8, "short row in PARASAIL_CASES: {line:.80}");
            let (s1, s2) = (f[0].as_bytes(), f[1].as_bytes());
            let want: i32 = f[3].parse().expect("score column");
            let sc = Scoring {
                match_score: f[4].parse().expect("match column"),
                mismatch: f[5].parse().expect("mismatch column"),
                open: f[6].parse().expect("open column"),
                ext: f[7].parse().expect("ext column"),
            };
            if s1.is_empty() || s2.is_empty() {
                continue;
            }
            n += 1;
            match semiglobal_ops(s1, s2, sc, &mut ops) {
                None => refused += 1,
                Some(got) => {
                    // identity from the recorded cigar
                    let (mut eq, mut total) = (0u32, 0u32);
                    let mut num = 0u32;
                    for c in f[2].bytes() {
                        if c.is_ascii_digit() {
                            num = num * 10 + (c - b'0') as u32;
                        } else {
                            total += num;
                            if c == b'=' {
                                eq += num;
                            }
                            num = 0;
                        }
                    }
                    let idy = if total > 0 {
                        eq as f64 / total as f64
                    } else {
                        0.0
                    };
                    let bi = if idy >= 0.95 {
                        4
                    } else if idy >= 0.90 {
                        3
                    } else if idy >= 0.80 {
                        2
                    } else if idy >= 0.60 {
                        1
                    } else {
                        0
                    };
                    if got < want {
                        worse += 1;
                        worst = worst.min(got - want);
                        buckets[bi][0] += 1;
                        let d = (want - got).unsigned_abs() as u64;
                        worse_by_delta[if d <= 2 {
                            0
                        } else if d <= 10 {
                            1
                        } else if d <= 50 {
                            2
                        } else {
                            3
                        }] += 1;
                        // Dump the first few real failures with their properties,
                        // since eleven crafted cases all agree and clearly are not
                        // reaching whatever distinguishes these.
                        // Every failure at >=95% identity: these are the only
                        // ones that can plausibly flip a merge verdict, since
                        // `align_to_merge` needs a >=100bp shared region within
                        // the diversity budget. Report the delta against that
                        // budget rather than in isolation.
                        if bi == 4 {
                            eprintln!(
                                "  HI-ID FAIL: len1={} len2={} idy={:.3} want={want} got={got} \
                                 delta={} cigar={:.50}",
                                s1.len(),
                                s2.len(),
                                idy,
                                got - want,
                                f[2]
                            );
                        }
                        if worse == 0 {
                            eprintln!(
                                "  FAIL#{worse}: len1={} len2={} want={want} got={got} \
                                 scoring=({}/{}/{}/{}) cigar={:.60}",
                                s1.len(),
                                s2.len(),
                                sc.match_score,
                                sc.mismatch,
                                sc.open,
                                sc.ext,
                                f[2]
                            );
                        }
                    } else if got > want {
                        better += 1;
                        buckets[bi][2] += 1;
                    } else {
                        buckets[bi][1] += 1;
                    }
                }
            }
        }
        eprintln!(
            "simd oracle: {n} cases | {worse} scored worse (worst delta {worst}) | \
             {better} scored higher | {refused} refused (exact path handles those)"
        );
        for (bi, label) in ["<60%", "60-80%", "80-90%", "90-95%", ">=95%"]
            .iter()
            .enumerate()
        {
            let [w, e, b] = buckets[bi];
            let tot = w + e + b;
            if tot > 0 {
                eprintln!(
                    "  identity {label:>7}: {tot:6} cases | {w:6} worse ({:.1}%) | {e:6} equal | {b:6} better",
                    100.0 * w as f64 / tot as f64
                );
            }
        }
        eprintln!(
            "  |delta| among the worse: 1-2={} 3-10={} 11-50={} 51+={}",
            worse_by_delta[0], worse_by_delta[1], worse_by_delta[2], worse_by_delta[3]
        );
        assert_eq!(
            worse, 0,
            "scored worse than exact parasail on {worse} of {n} cases (worst delta \
             {worst}); a different equally-optimal path is fine, a worse score is not"
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ops_span(ops: &[CigarOp]) -> (usize, usize) {
        let mut q = 0;
        let mut r = 0;
        for &(l, op) in ops {
            match op {
                b'=' | b'X' | b'M' => {
                    q += l;
                    r += l;
                }
                b'I' => q += l,
                b'D' => r += l,
                _ => {}
            }
        }
        (q, r)
    }

    #[test]
    fn identical_sequences_align_end_to_end() {
        let s = b"ACGTACGTACGTTTGCAAGGCTTAACCGGATTCAGGTACGATCGATCGGCTAACGT".repeat(4);
        let mut ops = Vec::new();
        assert!(semiglobal_ops(&s, &s, Scoring::MERGE, &mut ops).is_some());
        let (q, r) = ops_span(&ops);
        assert_eq!((q, r), (s.len(), s.len()));
        assert!(
            ops.iter().all(|&(_, o)| o == b'=' || o == b'M'),
            "identical input should produce only matches, got {ops:?}"
        );
    }

    #[test]
    fn a_contained_sequence_aligns_inside_the_longer_one() {
        // The case `align_to_merge` exists for: free end gaps let the shorter
        // sequence sit inside the longer without paying for the overhang.
        let long = b"ACGTTGCAAGGCTTAACCGGATTCAGGTACGATCGATCGGCTAACGTACGTACGTACGT".repeat(8);
        assert!(long.len() > 260, "fixture must be long enough to slice");
        let short = long[40..200].to_vec();
        let mut ops = Vec::new();
        assert!(semiglobal_ops(&long, &short, Scoring::MERGE, &mut ops).is_some());
        let (q, r) = ops_span(&ops);
        // Ops are reported relative to the caller's (s1, s2) = (long, short).
        assert_eq!(q, long.len(), "s1 fully spanned");
        assert_eq!(r, short.len(), "s2 fully spanned");
    }

    #[test]
    fn the_orientation_is_symmetric() {
        // Swapping the arguments must transpose I and D, not change the alignment.
        let long = b"ACGTTGCAAGGCTTAACCGGATTCAGGTACGATCGATCGGCTAACGT".repeat(8);
        let short = long[30..150].to_vec();
        let (mut a, mut b) = (Vec::new(), Vec::new());
        assert!(semiglobal_ops(&long, &short, Scoring::MERGE, &mut a).is_some());
        assert!(semiglobal_ops(&short, &long, Scoring::MERGE, &mut b).is_some());
        let (qa, ra) = ops_span(&a);
        let (qb, rb) = ops_span(&b);
        assert_eq!((qa, ra), (long.len(), short.len()));
        assert_eq!(
            (qb, rb),
            (short.len(), long.len()),
            "reversing the arguments must reverse which side each op consumes"
        );
    }

    #[test]
    fn which_ends_are_free() {
        // Establish the objective empirically instead of reasoning from the flag
        // names. Query = shorter. If the *reference's* overhangs are free, adding
        // junk to the long side must not change the score.
        let core = b"ACGTTGCAAGGCTTAACCGGATTCAGGTACGATCGATCGGCTAACGT".repeat(6);
        let mut ops = Vec::new();
        let base = semiglobal_ops(&core, &core, Scoring::MERGE, &mut ops).unwrap();

        let mut ref_pad = b"TTTTTTTTTTTTTTTTTTTT".to_vec();
        ref_pad.extend_from_slice(&core);
        ref_pad.extend_from_slice(b"GGGGGGGGGGGGGGGGGGGG");
        let padded_ref = semiglobal_ops(&ref_pad, &core, Scoring::MERGE, &mut ops).unwrap();

        let mut q_pad = b"TTTTTTTTTTTTTTTTTTTT".to_vec();
        q_pad.extend_from_slice(&core);
        let padded_query = semiglobal_ops(&core, &q_pad, Scoring::MERGE, &mut ops).unwrap();

        eprintln!("  same/same          score {base}");
        eprintln!("  long side padded    score {padded_ref}  (free if == {base})");
        eprintln!("  short side padded   score {padded_query}  (free if == {base})");
        assert_eq!(padded_ref, base, "the long side's overhangs should be free");
    }

    /// Where exactly does block-aligner disagree with exact parasail?
    ///
    /// The 54 884-case oracle says 20 000 score worse and 8 812 higher, uniformly
    /// across length buckets --- so not saturation and not block sizing. These are
    /// hand-built cases with a known answer, one property per case, to localise it.
    #[test]
    fn diagnose_against_exact_parasail() {
        use crate::parasail::{self, TieBreak};
        // The recorded cases use mismatch **-8** (`Scoring::BUBBLE`), not -2.
        // With mismatch that steep relative to gap open 12 / ext 1, the optimal
        // alignment prefers long indels over mismatches --- and long indels are
        // exactly where an adaptive *band* fails. The first version of this
        // diagnostic used MERGE (-2) and so never probed that regime.
        let sc = Scoring::BUBBLE;
        let base = b"ACGTTGCAAGGCTTAACCGGATTCAGGTACGATCGATCGGCTAACGTACGT".to_vec();

        let mut cases: Vec<(&str, Vec<u8>, Vec<u8>)> = Vec::new();
        cases.push(("identical", base.clone(), base.clone()));
        // one substitution in the middle
        let mut sub = base.clone();
        sub[25] = if sub[25] == b'A' { b'C' } else { b'A' };
        cases.push(("1 substitution", base.clone(), sub));
        // a 1-base deletion in the middle
        let mut del1 = base.clone();
        del1.remove(25);
        cases.push(("1-base gap", base.clone(), del1));
        // a 5-base deletion in the middle
        let mut del5 = base.clone();
        for _ in 0..5 {
            del5.remove(25);
        }
        cases.push(("5-base gap", base.clone(), del5));
        // a leading overhang (free in semi-global)
        let mut lead = b"TTTTTTTTTT".to_vec();
        lead.extend_from_slice(&base);
        cases.push(("leading overhang", lead, base.clone()));
        // a trailing overhang (free in semi-global)
        let mut trail = base.clone();
        trail.extend_from_slice(b"TTTTTTTTTT");
        cases.push(("trailing overhang", trail, base.clone()));
        // a non-ACGT byte: parasail's catch-all row scores it 0 against anything
        let mut ambig = base.clone();
        ambig[25] = b'N';
        cases.push(("one N", base.clone(), ambig));
        // Long internal indels: cheap under mismatch -8, and the case a band misses.
        for glen in [20usize, 60, 150] {
            // base is 51bp; repeat it to get room for a 150-base excision.
            let long_del: Vec<u8> = base.iter().cycle().take(600).copied().collect();
            let mut with_gap = long_del.clone();
            with_gap.drain(200..200 + glen);
            cases.push((
                match glen {
                    20 => "20-base internal gap",
                    60 => "60-base internal gap",
                    _ => "150-base internal gap",
                },
                long_del,
                with_gap,
            ));
        }
        // Many scattered mismatches: prefers gaps when mismatch is -8.
        let mut peppered = base.clone();
        for i in (5..base.len()).step_by(7) {
            peppered[i] = if peppered[i] == b'A' { b'C' } else { b'A' };
        }
        cases.push(("mismatch every 7bp", base.clone(), peppered));

        let mut ops = Vec::new();
        eprintln!(
            "  {:<20} {:>8} {:>8}   verdict",
            "case", "parasail", "block"
        );
        for (name, a, b) in cases {
            let want = parasail::semiglobal_with(&a, &b, sc, TieBreak::PARASAIL).score;
            let got = semiglobal_ops(&a, &b, sc, &mut ops);
            let g = got
                .map(|v| v.to_string())
                .unwrap_or_else(|| "refused".into());
            let verdict = match got {
                Some(v) if v == want => "same",
                Some(v) if v < want => "BLOCK LOWER",
                Some(_) => "BLOCK HIGHER",
                None => "refused",
            };
            eprintln!("  {name:<20} {want:>8} {g:>8}   {verdict}");
        }
    }

    #[test]
    fn empty_input_is_refused_rather_than_guessed() {
        let mut ops = Vec::new();
        assert!(semiglobal_ops(b"", b"ACGT", Scoring::MERGE, &mut ops).is_none());
        assert!(semiglobal_ops(b"ACGT", b"", Scoring::MERGE, &mut ops).is_none());
        assert!(ops.is_empty());
    }
}
