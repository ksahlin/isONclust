//! CARRIED ACROSS from the isONform port, which carried it from isONcorrect.
//! The algorithm below is unchanged; only this header is new.
//!
//! It is taken from **isONform**, not isONcorrect, deliberately: isONform's
//! finding 25 fixed two real defects in isONcorrect's copy after verifying it
//! against 56 549 recorded parasail calls -- the semi-global scan admitted an
//! "align nothing" end cell that parasail excludes, and broke the remaining
//! end-cell ties by a rule parasail does not use. isONcorrect's copy is 63 KB
//! and predates both fixes; this is the 79 KB one. Taking the older copy would
//! have been a silent regression.
//!
//! isONform's finding 24 is also in here and matters directly for isONclust: a
//! non-ACGT character scores **0** against itself, not as a match, because
//! `parasail.matrix_create("ACGT", ...)` only defines the four bases. isONclust
//! builds its matrix the same way, so an `N` in a read is scored 0 here too.
//! See PORTING.md's note on 2-bit encoding.
//!
//! WHAT ISONCLUST USES IT FOR, and how the call differs from isONcorrect's:
//!
//! ```text
//! isONcorrect  match  4, mismatch -8, open 12, ext 1   (the overcorrection guard)
//! isONclust    match  2, mismatch -2, open 2..5, ext 1 (get_best_cluster_block_align)
//! ```
//!
//! The opening penalty is chosen per comparison from the summed error rates of
//! the two reads, so all four values occur. Same function, same tie-breaking,
//! different `Scoring` -- which is why this module is parameterised rather than
//! hard-coded.
//!
//! Semi-global affine alignment matching `parasail.sg_trace_scan_16`.
//!
//! This is the alignment behind the structural-overcorrection guard: the
//! corrected read is aligned back against the original, and [`crate::guard`]
//! reverts any indel run longer than 10 columns on the assumption it is real
//! structure (an exon difference) rather than an error.
//!
//! `correct_read` calls it as
//!
//! ```text
//! parasail_alignment(seq, corr, match_score=4, mismatch_penalty=-8,
//!                    opening_penalty=12, gap_ext=1)
//! ```
//!
//! — note 12, not the function's own default of 24.
//!
//! # What "sg" means, measured
//!
//! parasail's plain `sg` leaves gaps at **both ends of both sequences** free.
//! Verified against the library rather than read off documentation:
//!
//! ```text
//! sg("ACGTACGT", "TTTACGTACGTTTT") -> 3D8=3D, score 32 == 8 matches * 4
//! sg("TTTACGTACGTTTT", "ACGTACGT") -> 3I8=3I, score 32
//! ```
//!
//! Both end gaps contribute nothing to the score, and the CIGAR still spans
//! both sequences in full. A gap of length `L` inside the alignment costs
//! `open + (L - 1) * ext`, so a single-base gap costs `open` — also measured:
//! `sg("ACGTACGT", "ACGACGT")` scores 16, which is the free-end-gap alignment
//! tying exactly with `7 * 4 - 12` for the one-base-deletion alignment.
//!
//! That tie is the whole problem. Free end gaps plus affine costs leave many
//! alignments at the optimal score, and parasail returns one of them. Which one
//! changes where indel runs fall, which changes what the guard reverts, which
//! changes the output — so the choice has to be reproduced, not merely matched
//! on score.
//!
//! # CIGAR convention
//!
//! `I` consumes `s1` (the query), `D` consumes `s2` (the database), matching
//! what `cigar_to_seq(cigar, s1, s2)` then expects — the same convention as
//! edlib's, so [`crate::align::cigar_to_seq`] expands these unchanged.
//!
//! `#![allow(dead_code)]` below is deliberate: this file is carried across
//! VERBATIM from the isONform port and is verified there. isONclust uses only
//! part of its surface today, and trimming the rest would mean editing code
//! this port has not verified, to silence a lint about code that is correct.
//! The unused parts are what the remaining stages will call.
#![allow(dead_code)]

use crate::align::CigarOp;

/// Scoring, in parasail's sign convention: `match_score` positive, `mismatch`
/// negative, and both gap penalties positive numbers that get subtracted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Scoring {
    pub match_score: i32,
    pub mismatch: i32,
    pub open: i32,
    pub ext: i32,
}

impl Scoring {
    /// What `align_bubble_nodes` passes when deciding whether a bubble is
    /// poppable (`SimplifyGraph.py:656`): `match=2, mismatch=-8, open=12, ext=1`.
    ///
    /// Note the mismatch penalty overrides `parasail_alignment`'s own default of
    /// −2. isONcorrect's equivalent was `match=4, mismatch=-8`; the constant it
    /// used (`GUARD`, for its overcorrection guard) is gone, because neither the
    /// guard nor that scoring exists here.
    pub const BUBBLE: Scoring = Scoring {
        match_score: 2,
        mismatch: -8,
        open: 12,
        ext: 1,
    };

    /// What `IsoformGeneration.align_to_merge` passes (`IsoformGeneration.py:379`),
    /// which is also `parasail_alignment`'s own default: `match=2, mismatch=-2,
    /// open=12, ext=1`.
    ///
    /// The two differ only in the mismatch penalty, and that difference is the
    /// whole distinction between "are these two bubble paths the same" and "are
    /// these two isoforms the same" --- worth not conflating.
    pub const MERGE: Scoring = Scoring {
        match_score: 2,
        mismatch: -2,
        open: 12,
        ext: 1,
    };
}

/// Which way a tie is broken. Every field was determined by sweeping against
/// recorded parasail output; see the `sweep` test.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TieBreak {
    /// Order to try the three predecessors of an `H` cell.
    pub h_order: [Step; 3],
    /// On a tie inside a gap chain, prefer opening a new gap over extending.
    pub prefer_open: bool,
    /// When several cells on the last row/column share the best score, scan the
    /// last column before the last row.
    pub column_first: bool,
    /// Keep the last equally-good end cell rather than the first.
    pub last_max: bool,
}

/// A predecessor of an `H` cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Step {
    /// Diagonal: a match or mismatch.
    Diagonal,
    /// From the gap-in-`s1` chain, emitting `D`.
    Delete,
    /// From the gap-in-`s2` chain, emitting `I`.
    Insert,
}

impl TieBreak {
    /// The ordering that reproduces `sg_trace_scan_16`, measured rather than
    /// guessed — see the `sweep` test for the numbers.
    ///
    /// `column_first` and `last_max` only matter when several cells on the last
    /// row or column tie for the best score. On isONcorrect's corpus that never
    /// happened — a read against its own correction always prefers a full-length
    /// alignment — so both were **unpinned by evidence**, with the sweep reporting
    /// all four combinations as equally perfect and a note to re-run it if a
    /// corpus ever reached such a tie.
    ///
    /// **isONform's does, constantly**, because it aligns two bubble-path
    /// consensuses rather than a read against itself. Swept over 54 884 recorded
    /// isONform calls the winner is now unique and strict:
    ///
    /// ```text
    ///   54884 / 54884  [Diagonal, Insert, Delete] open=false col_first=false last_max=false
    ///   54664 / 54884  [Diagonal, Insert, Delete] open=false col_first=true  last_max=false
    ///   52891 / 54884  [Diagonal, Insert, Delete] open=false col_first=false last_max=true
    /// ```
    ///
    /// So all four fields are pinned by evidence now. Note this only became a
    /// clean sweep after `end_cell` stopped admitting the corner into the row and
    /// column ranges (finding 25); before that the best any setting reached was
    /// 54 772.
    pub const PARASAIL: TieBreak = TieBreak {
        h_order: [Step::Diagonal, Step::Insert, Step::Delete],
        prefer_open: false,
        column_first: false,
        last_max: false,
    };

    /// Every combination, for the sweep that picks the real one.
    #[cfg(test)]
    fn all() -> Vec<TieBreak> {
        use Step::*;
        let orders = [
            [Diagonal, Delete, Insert],
            [Diagonal, Insert, Delete],
            [Delete, Diagonal, Insert],
            [Delete, Insert, Diagonal],
            [Insert, Diagonal, Delete],
            [Insert, Delete, Diagonal],
        ];
        let mut out = Vec::new();
        for h_order in orders {
            for prefer_open in [true, false] {
                for column_first in [true, false] {
                    for last_max in [true, false] {
                        out.push(TieBreak {
                            h_order,
                            prefer_open,
                            column_first,
                            last_max,
                        });
                    }
                }
            }
        }
        out
    }
}

/// Result of a semi-global alignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    pub score: i32,
    pub cigar: String,
    pub ops: Vec<CigarOp>,
}

const NEG: i32 = i32::MIN / 4;

/// parasail's alphabet classes for `parasail_matrix_create("ACGT", …)`: one index
/// per base plus a single catch-all for every other byte. Case-insensitive,
/// because parasail's own mapping is (`matrix_create` accepts `a` for `A`).
const CATCHALL: u8 = 4;

const CLASS: [u8; 256] = {
    let mut t = [CATCHALL; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t
};

/// The substitution score for two characters, as
/// `parasail_matrix_create("ACGT", match, mismatch)` defines it.
///
/// **Not byte equality**, and the gap between the two is reachable rather than
/// theoretical. `matrix_create` builds a 5×5 matrix — one row and column per
/// base, plus one catch-all — and fills the catch-all row and column with
/// **zero**. So any byte outside `ACGT` scores 0 against everything, *including
/// an identical copy of itself*: real parasail scores `X` against `X` as 0, not
/// as a match. Measured, not inferred; see the tests below.
///
/// That is a live code path, not a curiosity. `generate_consensus_path` returns
/// `"X" * max_len` whenever every path span in a bubble is shorter than `k`, so
/// two such placeholders get compared, and scoring X-against-X as a match is
/// exactly the case byte equality gets wrong. It made the port pop bubbles the
/// reference refuses --- `PORTING.md`, finding 24.
///
/// The CIGAR *operation* still comes from byte equality, which is not an
/// inconsistency: parasail emits `=` for two identical characters even when the
/// matrix scores them 0, so `X` against `X` is a zero-scoring `=`.
#[inline(always)]
fn subst(sc: Scoring, a: u8, b: u8) -> i32 {
    let (ca, cb) = (CLASS[a as usize], CLASS[b as usize]);
    if ca == CATCHALL || cb == CATCHALL {
        0
    } else if ca == cb {
        sc.match_score
    } else {
        sc.mismatch
    }
}

/// Align `s1` against `s2`, reproducing `parasail.sg_trace_scan_16`.
///
/// `O(n*m)` time and memory. The guard runs this once per read against a
/// sequence of nearly the same length, so the table is the read length squared —
/// the one place in the port where that is a real cost, and the reason
/// *Deferred improvements* notes banding it.
pub fn semiglobal(s1: &[u8], s2: &[u8], sc: Scoring) -> Alignment {
    semiglobal_with(s1, s2, sc, TieBreak::PARASAIL)
}

/// The four traceback bits per cell, plus the scores the end-cell scan needs.
///
/// # Why the decisions are stored rather than the matrices
///
/// The obvious implementation keeps all three `i32` matrices and re-derives each
/// traceback step by comparing them. That is 12 bytes per cell — 12 MB for a
/// 1 000 bp read against its own correction, once per read — and it was **65% of
/// the port's runtime** once the edit-distance hot spot was fixed.
///
/// Two rows of `H`/`E`/`F` are enough for the forward pass. What the traceback
/// actually needs is not the scores but the **choices**, and there are only four
/// bits of those per cell:
///
/// * which of the three predecessors won the `H` cell (2 bits), decided by
///   [`TieBreak::h_order`];
/// * whether a `D` chain leaves for `H` here or extends (1 bit);
/// * the same for an `I` chain (1 bit).
///
/// One byte per cell that is `n * m` — 1 MB instead of 12 MB — and it is the
/// *same* computation, so the recorded CIGARs must not move.
///
/// Four bits would halve it again, and was tried: the read-modify-write and
/// shifting cost about 11% of runtime while the DP is already streaming
/// row-by-row, so the cache never benefited. One byte per cell is the better
/// trade — measured, not assumed.
struct Table {
    m: usize,
    packed: Vec<u8>,
    /// `H[i][m]` for every row, for the last-column scan.
    last_col: Vec<i32>,
    /// `H[n][j]` for every column, for the last-row scan.
    last_row: Vec<i32>,
}

impl Table {
    #[inline]
    fn get(&self, i: usize, j: usize) -> u8 {
        self.packed[(i - 1) * self.m + (j - 1)]
    }

    #[inline]
    fn h_choice(&self, i: usize, j: usize) -> Step {
        match self.get(i, j) & 0b11 {
            0 => Step::Diagonal,
            1 => Step::Delete,
            _ => Step::Insert,
        }
    }

    #[inline]
    fn delete_leaves(&self, i: usize, j: usize) -> bool {
        self.get(i, j) & 0b100 != 0
    }

    #[inline]
    fn insert_leaves(&self, i: usize, j: usize) -> bool {
        self.get(i, j) & 0b1000 != 0
    }
}

/// Row bounds of the band, clamped to the matrix.
#[inline]
fn band(i: usize, half: usize, m: usize) -> (usize, usize) {
    (i.saturating_sub(half).max(1), (i + half).min(m))
}

/// Forward pass: fill the decision table row by row, keeping only two rows of
/// scores alive.
///
/// `half` restricts each row to `|i - j| <= half`; `None` computes every cell.
/// Cells outside the band keep `NEG`, so they can never win the end-cell scan.
fn forward(
    s1: &[u8],
    s2: &[u8],
    sc: Scoring,
    tb: TieBreak,
    half: Option<usize>,
    global: bool,
) -> Table {
    let (n, m) = (s1.len(), s2.len());

    // Free end gaps: an unaligned prefix of either sequence costs nothing, so
    // row 0 and column 0 start at zero — including the gap chains, which is
    // what keeps a *leading* gap free rather than charging `open`.
    let mut h_prev = vec![0i32; m + 1];
    let mut e_prev = vec![0i32; m + 1];
    let mut f_prev = vec![NEG; m + 1];
    f_prev[0] = 0;
    if global {
        // Row 0 of a *global* alignment is a gap run, charged. Semi-global leaves
        // it at zero, which is what makes a leading gap free.
        for j in 1..=m {
            h_prev[j] = -(sc.open + sc.ext * (j as i32 - 1));
            e_prev[j] = h_prev[j];
        }
    }

    let mut h_cur = vec![NEG; m + 1];
    let mut e_cur = vec![NEG; m + 1];
    let mut f_cur = vec![NEG; m + 1];

    let mut table = Table {
        m,
        packed: vec![0u8; n * m],
        last_col: vec![NEG; n + 1],
        last_row: vec![NEG; m + 1],
    };
    table.last_col[0] = 0;

    // Mask of achieving predecessors -> the winner under `tb.h_order`. Bit 0 is
    // the diagonal, 1 is delete, 2 is insert. Mask 0 cannot happen, since the
    // maximum is achieved by something.
    let choice_lut: [u8; 8] = std::array::from_fn(|mask| {
        for step in tb.h_order {
            let (bit, code) = match step {
                Step::Diagonal => (0, 0u8),
                Step::Delete => (1, 1),
                Step::Insert => (2, 2),
            };
            if mask >> bit & 1 == 1 {
                return code;
            }
        }
        0
    });

    for i in 1..=n {
        let (lo, hi) = match half {
            Some(h) => band(i, h, m),
            None => (1, m),
        };
        if half.is_some() {
            // Outside the band a cell is unreachable, not zero.
            h_cur.fill(NEG);
            e_cur.fill(NEG);
            f_cur.fill(NEG);
        }
        if lo == 1 {
            h_cur[0] = if global {
                -(sc.open + sc.ext * (i as i32 - 1))
            } else {
                0
            };
            e_cur[0] = NEG;
            f_cur[0] = if global { h_cur[0] } else { 0 };
        }
        let c1 = s1[i - 1];
        // One row of decisions, sliced out once: the index arithmetic and bounds
        // check would otherwise be repeated for every cell.
        let row = &mut table.packed[(i - 1) * m..i * m];

        for j in lo..=hi {
            let open_d = h_cur[j - 1] - sc.open;
            let ext_d = e_cur[j - 1] - sc.ext;
            e_cur[j] = open_d.max(ext_d);
            let d_leaves = if tb.prefer_open {
                open_d >= ext_d
            } else {
                open_d > ext_d
            };

            let open_i = h_prev[j] - sc.open;
            let ext_i = f_prev[j] - sc.ext;
            f_cur[j] = open_i.max(ext_i);
            let i_leaves = if tb.prefer_open {
                open_i >= ext_i
            } else {
                open_i > ext_i
            };

            let diag = h_prev[j - 1] + subst(sc, c1, s2[j - 1]);
            let best = diag.max(e_cur[j]).max(f_cur[j]);
            h_cur[j] = best;

            // Which predecessors achieve the cell, as a 3-bit mask resolved
            // through the table above. Walking `tb.h_order` with an early break
            // instead is a data-dependent branch chain over a runtime array,
            // executed once per cell, and the predictor cannot learn it.
            let mask = usize::from(diag == best)
                | usize::from(e_cur[j] == best) << 1
                | usize::from(f_cur[j] == best) << 2;
            row[j - 1] = choice_lut[mask] | (u8::from(d_leaves) << 2) | (u8::from(i_leaves) << 3);
        }

        // The last column is only reachable when the band reaches it.
        table.last_col[i] = if hi == m { h_cur[m] } else { NEG };
        std::mem::swap(&mut h_prev, &mut h_cur);
        std::mem::swap(&mut e_prev, &mut e_cur);
        std::mem::swap(&mut f_prev, &mut f_cur);
    }

    // After the last swap `h_prev` holds row n --- or row 0 when the loop never
    // ran, which is all zeros and correct.
    table.last_row.copy_from_slice(&h_prev);
    table
}

/// [`forward`], evaluated along anti-diagonals instead of rows.
///
/// # Why this order, and why it is exact
///
/// The row-wise loop has a serial dependency: `e_cur[j]` reads `h_cur[j - 1]` and
/// `e_cur[j - 1]`, so a row cannot be filled in parallel. That is the dependency
/// Farrar's striped SIMD works around with a lazy-F correction pass, and getting
/// the *per-cell tie-break byte* bit-exact through that correction is delicate.
///
/// On an anti-diagonal `d = i + j` there is no such dependency. Every input a
/// cell needs sits on `d - 1` or `d - 2`:
///
/// | needed | cell | diagonal |
/// | --- | --- | --- |
/// | `H(i, j-1)`, `E(i, j-1)` | left | `d - 1` |
/// | `H(i-1, j)`, `F(i-1, j)` | above | `d - 1` |
/// | `H(i-1, j-1)` | diagonal | `d - 2` |
///
/// So the cells of one diagonal are mutually independent and the recurrence is
/// unchanged --- same arithmetic, same comparisons, same tie-break, evaluated in
/// a different but equally valid order. That makes this a stepping stone rather
/// than a rewrite: it is verified against [`forward`] by
/// `wavefront_matches_the_row_order_exactly` and by the `parasail::oracle`, and
/// the SIMD kernel then only has to match *this*.
///
/// Indexing everything by row `i` is what keeps it readable: on diagonal `d`,
/// `prev1[i]` is the cell to the left and `prev1[i - 1]` the cell above.
///
/// Scalar, this is *slower* than [`forward`] --- strided access, and a diagonal's
/// length varies. Its value is that it is exact and parallel.
fn forward_wavefront(s1: &[u8], s2: &[u8], sc: Scoring, tb: TieBreak) -> Table {
    let (n, m) = (s1.len(), s2.len());
    let mut table = Table {
        m,
        packed: vec![0u8; n * m],
        last_col: vec![NEG; n + 1],
        last_row: vec![NEG; m + 1],
    };
    table.last_col[0] = 0;
    // The degenerate shapes, matching `forward` exactly. With `m == 0` its row
    // loop still runs and sets `h_cur[0] = 0` before finding an empty inner
    // range, so every last-column entry is 0, not NEG. With `n == 0` the row
    // loop never runs and the last row stays row 0, which is all zeros.
    if m == 0 {
        table.last_col.fill(0);
        table.last_row[0] = 0;
        return table;
    }
    if n == 0 {
        table.last_row.fill(0);
        return table;
    }

    let choice_lut = choice_lut(tb);

    // Indexed by row, 0..=n. Diagonal d holds rows max(1, d-m)..=min(n, d-1).
    // The row-0 and column-0 boundaries are seeded exactly as `forward` does:
    // H is 0 on both, E is 0 on row 0 and NEG on column 0, F is NEG on row 0 and
    // 0 on column 0.
    let mut p2_h = vec![NEG; n + 1];
    let mut p1_h = vec![NEG; n + 1];
    let mut p1_e = vec![NEG; n + 1];
    let mut p1_f = vec![NEG; n + 1];
    let mut cur_h = vec![NEG; n + 1];
    let mut cur_e = vec![NEG; n + 1];
    let mut cur_f = vec![NEG; n + 1];
    let mut scratch = vec![0i32; n + 1];
    // `s2` reversed once, so the diagonal's backwards walk over it is a
    // contiguous forward read.
    let s2_rev: Vec<u8> = s2.iter().rev().copied().collect();

    // d = 0: only (0, 0).
    p2_h[0] = 0;
    // d = 1: (0, 1) and (1, 0), both boundary cells.
    p1_h[0] = 0; // (0, 1)
    p1_e[0] = 0; // E(0, j) = 0
    p1_f[0] = NEG; // F(0, j) = NEG for j >= 1
    if n >= 1 {
        p1_h[1] = 0; // (1, 0)
        p1_e[1] = NEG; // E(i, 0) = NEG
        p1_f[1] = 0; // F(i, 0) = 0
    }

    for d in 2..=(n + m) {
        let lo = d.saturating_sub(m).max(1);
        let hi = n.min(d - 1);
        // The boundary cells of this diagonal, if it touches an edge.
        if d <= m {
            cur_h[0] = 0;
            cur_e[0] = 0;
            cur_f[0] = NEG;
        }
        if d <= n {
            cur_h[d] = 0;
            cur_e[d] = NEG;
            cur_f[d] = 0;
        }
        // Split into phases over contiguous slices so the arithmetic is
        // vectorisable: `e`, `f` and `h` are elementwise over runs of `i`, and
        // only the packed-byte store is scattered (its index strides by `m - 1`).
        let width = hi + 1 - lo;
        let e = &mut cur_e[lo..=hi];
        let f = &mut cur_f[lo..=hi];

        // E(i, j) from the cell to the left: both inputs at row `i` on d-1.
        for (k, ev) in e.iter_mut().enumerate() {
            let i = lo + k;
            *ev = (p1_h[i] - sc.open).max(p1_e[i] - sc.ext);
        }
        // F(i, j) from the cell above: both inputs at row `i - 1` on d-1.
        for (k, fv) in f.iter_mut().enumerate() {
            let i = lo + k;
            *fv = (p1_h[i - 1] - sc.open).max(p1_f[i - 1] - sc.ext);
        }
        // The diagonal score. `j = d - i` runs backwards as `i` runs forwards, so
        // `s2` is read in reverse --- taken from a reversed copy made once, which
        // keeps this a contiguous load rather than a gather.
        let s2r_from = m - (d - lo);
        for k in 0..width {
            let i = lo + k;
            scratch[k] = p2_h[i - 1] + subst(sc, s1[i - 1], s2_rev[s2r_from + k]);
        }
        for k in 0..width {
            let i = lo + k;
            cur_h[i] = scratch[k].max(e[k]).max(f[k]);
        }

        // Tie-break bytes. Scattered store, so scalar --- but the comparisons are
        // over values already in cache from the phases above.
        for k in 0..width {
            let i = lo + k;
            let j = d - i;
            let best = cur_h[i];
            let (ev, fv) = (e[k], f[k]);
            let mask = usize::from(scratch[k] == best)
                | usize::from(ev == best) << 1
                | usize::from(fv == best) << 2;
            let d_leaves = if tb.prefer_open {
                p1_h[i] - sc.open >= p1_e[i] - sc.ext
            } else {
                p1_h[i] - sc.open > p1_e[i] - sc.ext
            };
            let i_leaves = if tb.prefer_open {
                p1_h[i - 1] - sc.open >= p1_f[i - 1] - sc.ext
            } else {
                p1_h[i - 1] - sc.open > p1_f[i - 1] - sc.ext
            };
            table.packed[(i - 1) * m + (j - 1)] =
                choice_lut[mask] | (u8::from(d_leaves) << 2) | (u8::from(i_leaves) << 3);
            if j == m {
                table.last_col[i] = best;
            }
            if i == n {
                table.last_row[j] = best;
            }
        }
        std::mem::swap(&mut p2_h, &mut p1_h);
        std::mem::swap(&mut p1_h, &mut cur_h);
        std::mem::swap(&mut p1_e, &mut cur_e);
        std::mem::swap(&mut p1_f, &mut cur_f);
        cur_h.fill(NEG);
        cur_e.fill(NEG);
        cur_f.fill(NEG);
    }
    table.last_row[0] = 0;
    table
}

/// The mask-to-winner table, shared by both fill orders.
fn choice_lut(tb: TieBreak) -> [u8; 8] {
    std::array::from_fn(|mask| {
        for step in tb.h_order {
            let (bit, code) = match step {
                Step::Diagonal => (0, 0u8),
                Step::Delete => (1, 1),
                Step::Insert => (2, 2),
            };
            if mask >> bit & 1 == 1 {
                return code;
            }
        }
        0
    })
}

/// [`semiglobal`] with an explicit tie-break, for the sweep.
///
/// # Banding: measured, and rejected
///
/// The full DP is `O(n*m)` and is the largest single cost in the port. The two
/// sequences are a read and its own correction, so the optimal path hugs the
/// diagonal — over 2 408 recorded guard alignments the largest excursion is 63
/// cells, and a half-band of 64 reproduces **every** parasail CIGAR exactly.
/// The `banding` tests below are that measurement.
///
/// It is still not used, for two reasons that only became clear once both were
/// measured:
///
/// * **the cheap band is not provably correct.** There is a real bound — every
///   column is diagonal or a gap, so `2c + g = n + m`, and reaching deviation
///   `W` needs `g >= W`, hence any such path scores at most `2(n+m) - 2W`. But
///   the scores here sit far below the all-match maximum (median slack 753), so
///   the band that *proves* optimality is `~slack/2 ≈ 380`, not 64. And even
///   that proves only the score: a cell near the band edge can have an
///   under-estimated `E`/`F`, which could flip a tie-break and change the CIGAR;
/// * **the provable band is worth almost nothing.** Measured on a 200-read gene
///   cluster: full DP 3.4 s, provable band 3.1 s, unproven ±64 band 2.1 s. The
///   rigorous version buys 9%.
///
/// So the choice was 9% with a proof, or 35% on evidence that would have to be
/// re-established for every new kind of input. For a tool whose entire contract
/// is byte-identical output, and whose failure mode here is silent, neither is a
/// good trade. The remaining win is vectorising this loop, which is exact
/// because it does not change the recurrence.
/// Global (Needleman-Wunsch) alignment: end gaps are **charged**.
///
/// # Why this exists
///
/// `align_bubble_nodes` aligns the two paths of a bubble, and those paths are
/// delimited by the bubble's *shared* start and end nodes --- `con =
/// all_reads[q_id][1][pos1 : pos2 + k_size]`, with `pos1`/`pos2` the read's
/// positions at those anchors. Sequences with common endpoints are a global
/// alignment problem, but the reference calls `parasail.sg_trace_scan_16`
/// (semi-global), which makes both end gaps free.
///
/// That matters for the decision it feeds. `mergeable_start`/`mergeable_end`
/// measure how much dangles outside the first and last significant match, bounded
/// by `delta_iso_len_5`/`_3`. If the aligner may shed the ends for nothing, it
/// will prefer alignments that look well-anchored by trimming exactly the
/// end disagreement those thresholds exist to catch --- biasing toward popping
/// bubbles whose paths differ at their boundaries. The reference's separate
/// `(longer_len - shorter_len) / longer_len < delta` gate at `delta = 0.20` is
/// arguably compensation for that.
///
/// Selected by `ISONFORM_BUBBLE_GLOBAL=1`. Off by default: it changes which
/// bubbles pop, so it is a method change to be measured, not assumed.
pub fn global(s1: &[u8], s2: &[u8], sc: Scoring) -> Alignment {
    global_with(s1, s2, sc, TieBreak::PARASAIL)
}

/// [`global`] with an explicit tie-break.
///
/// The only difference from [`semiglobal_with`] is the boundary condition: row 0
/// and column 0 are gap-penalised rather than free, and the traceback starts at
/// `(n, m)` rather than the best cell on the last row or column.
pub fn global_with(s1: &[u8], s2: &[u8], sc: Scoring, tb: TieBreak) -> Alignment {
    let (n, m) = (s1.len(), s2.len());
    if n == 0 || m == 0 {
        // One side empty: the whole of the other is a single gap run.
        let (len, op) = if n == 0 { (m, b'D') } else { (n, b'I') };
        let score = if len == 0 {
            0
        } else {
            -(sc.open + sc.ext * (len as i32 - 1))
        };
        return Alignment {
            score,
            cigar: if len == 0 {
                String::new()
            } else {
                format!("{len}{}", op as char)
            },
            ops: if len == 0 { vec![] } else { vec![(len, op)] },
        };
    }
    let table = forward(s1, s2, sc, tb, None, true);
    let mut aln = traceback_from(s1, s2, &table, n, m);
    aln.score = table.last_row[m];
    aln
}

pub fn semiglobal_with(s1: &[u8], s2: &[u8], sc: Scoring, tb: TieBreak) -> Alignment {
    // Every cell, always: banding was implemented, measured and rejected above.
    //
    // Which fill order runs is a one-time read, not a per-call branch. The
    // wavefront is exact --- verified cell-for-cell against the row order and on
    // 54 884 real-library cases --- and exists as the base the SIMD kernel is
    // built and checked against. It is *slower* scalar, so the row order stays
    // the default until the vectorised kernel lands.
    if wavefront_enabled() {
        traceback(s1, s2, tb, &forward_wavefront(s1, s2, sc, tb))
    } else {
        traceback(s1, s2, tb, &forward(s1, s2, sc, tb, None, false))
    }
}

/// `ISONFORM_WAVEFRONT=1` selects the anti-diagonal fill. Read once.
fn wavefront_enabled() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("ISONFORM_WAVEFRONT").is_some())
}

// ---------------------------------------------------------------------------
// Removed on the way across from isONcorrect
// ---------------------------------------------------------------------------
//
// Two blocks of research scaffolding were dropped here, both isONcorrect-only:
//
//  * `global_affine`, `global_check` and `global_check_report` with their
//    counters, which existed to answer whether exact affine *global* alignment
//    could replace this semi-global one. That question was about isONcorrect's
//    structural-overcorrection guard --- the comparison was whether
//    `guard::fix_correction` landed on the same corrected sequence --- and
//    isONform has no guard, so there is nothing to compare.
//  * `band_check` and `band_check_report`, which cross-checked a banded
//    alignment against the exact one on arbitrary real data. Banding was
//    measured and rejected in the isONcorrect port; the exact path is what
//    shipped, and this file keeps it.
//
// Both ran only when an `ISONCORRECT_*` environment variable was set, so
// dropping them changes nothing on any normal run. If banding is ever revisited
// here, that scaffolding is worth copying back rather than rewriting.

impl Table {
    /// The end cell the alignment starts its traceback from: the best-scoring cell
    /// on the last row or last column, since trailing gaps are free.
    ///
    /// Two things about parasail's choice here are not guessable, and both were
    /// measured rather than reasoned about. `result.end_query` / `result.end_ref`
    /// report the cell parasail picked, which is what made them measurable.
    ///
    /// **1. The scan starts at 1, not 0.** Cell `(n, 0)` is reached by consuming
    /// all of `s1` as a free leading gap and none of `s2`; `(0, m)` is its mirror.
    /// Both score 0, so including them lets "align nothing" win whenever every
    /// real alignment scores below zero. parasail does not do that --- it insists
    /// on at least one diagonal step:
    ///
    /// ```text
    /// m = parasail.matrix_create("ACGT", 2, -2)
    /// sg("A",    "C")    -> -2  "1X"
    /// sg("AAAA", "CCCC") -> -2  "3I1X3D"
    /// ```
    ///
    /// **2. The corner is excluded from both ranges and considered last.** That
    /// looks like a detail and decides real cases. All-A against all-C ties
    /// everywhere, so it isolates the rule; the cell parasail picks, 1-based:
    ///
    /// ```text
    ///   n\m     1      2      3      4      5      6
    ///     1   (1,1)  (1,1)  (1,1)  (1,1)  (1,1)  (1,1)
    ///     2   (1,1)  (2,1)  (2,1)  (2,1)  (2,1)  (2,1)
    ///     3   (1,1)  (3,1)  (3,1)  (3,1)  (3,1)  (3,1)
    ///     4   (1,1)  (4,1)  (4,1)  (4,1)  (4,1)  (4,1)
    /// ```
    ///
    /// The ties are always between `(n, 1)` on the last row and `(1, m)` on the
    /// last column. For `m >= 2` parasail takes the row one; for `m == 1` it takes
    /// `(1, 1)`, which is on the column. No row-before-column or
    /// column-before-row rule produces both — but excluding the corner does:
    /// with `m == 1` the row range `1..m` is **empty**, so the choice falls to the
    /// column and lands on its first maximum, `(1, 1)`.
    ///
    /// Verified exactly on **56 549** recorded real calls (54 884 Drosophila,
    /// 1 665 `sirv_small`): 0 score and 0 CIGAR mismatches. Before this, 112 and
    /// 35 CIGARs differed --- all equally *optimal*, which is why the score gate
    /// could not see them, and `parse_cigar_diversity` reads the CIGAR rather than
    /// the score. `PORTING.md` finding 25.
    ///
    /// Only gaps in row 0 / column 0 and gaps after the end cell are free; gaps in
    /// between are charged, which is what makes one mismatch cheaper than walking
    /// the border.
    fn end_cell(&self, tb: TieBreak) -> (i32, usize, usize) {
        let n = self.last_col.len() - 1;
        let m = self.m;
        // parasail rejects an empty input outright ("s1Len must be > 0"), so the
        // reference raises rather than returning anything here, and `decide`'s
        // length gate keeps it unreachable. Defined rather than left to panic on
        // an all-NEG scan.
        if n == 0 || m == 0 {
            return (0, n, m);
        }
        let mut best: (i32, usize, usize) = (NEG, n, m);
        // The corner is excluded from both ranges and considered last. With
        // `m == 1` that leaves the row scan empty, which is what makes parasail
        // pick `(1, 1)` there and `(n, 1)` for every `m >= 2` --- see the table in
        // the doc comment.
        let row = (1..m).map(|j| (self.last_row[j], n, j));
        let col = (1..n).map(|i| (self.last_col[i], i, m));
        let corner = (self.last_row[m], n, m);
        let consider = |cand: (i32, usize, usize), best: &mut (i32, usize, usize)| {
            if cand.0 > best.0 || (tb.last_max && cand.0 == best.0) {
                *best = cand;
            }
        };
        if tb.column_first {
            for c in col {
                consider(c, &mut best);
            }
            for r in row {
                consider(r, &mut best);
            }
        } else {
            for r in row {
                consider(r, &mut best);
            }
            for c in col {
                consider(c, &mut best);
            }
        }
        consider(corner, &mut best);
        best
    }
}

/// Walk the stored decisions back from the end cell.
fn traceback(s1: &[u8], s2: &[u8], tb: TieBreak, table: &Table) -> Alignment {
    let (score, i, j) = table.end_cell(tb);
    let mut aln = traceback_from(s1, s2, table, i, j);
    aln.score = score;
    aln
}

/// The traceback walk itself, from an explicit end cell.
fn traceback_from(s1: &[u8], s2: &[u8], table: &Table, i0: usize, j0: usize) -> Alignment {
    let (n, m) = (s1.len(), s2.len());
    let (mut i, mut j) = (i0, j0);

    let mut ops: Vec<u8> = Vec::with_capacity(n + m);
    // Free trailing gaps first: whatever is left of either sequence.
    ops.extend(std::iter::repeat_n(b'D', m - j));
    ops.extend(std::iter::repeat_n(b'I', n - i));

    // `Some(Step::Delete)` means we are inside a D chain and must decide
    // open-vs-extend before leaving it; `None` means we are on H.
    let mut chain: Option<Step> = None;
    while i > 0 && j > 0 {
        match chain {
            Some(Step::Delete) => {
                let leaves = table.delete_leaves(i, j);
                ops.push(b'D');
                j -= 1;
                chain = if leaves { None } else { Some(Step::Delete) };
            }
            Some(Step::Insert) => {
                let leaves = table.insert_leaves(i, j);
                ops.push(b'I');
                i -= 1;
                chain = if leaves { None } else { Some(Step::Insert) };
            }
            Some(Step::Diagonal) => unreachable!("the diagonal is not a chain"),
            None => match table.h_choice(i, j) {
                Step::Diagonal => {
                    ops.push(if s1[i - 1] == s2[j - 1] { b'=' } else { b'X' });
                    i -= 1;
                    j -= 1;
                }
                step => chain = Some(step),
            },
        }
    }
    // Free leading gaps.
    ops.extend(std::iter::repeat_n(b'D', j));
    ops.extend(std::iter::repeat_n(b'I', i));

    ops.reverse();
    let (cigar, runs) = encode(&ops);
    // The caller fills in the score; the walk does not know it.
    Alignment {
        score: 0,
        cigar,
        ops: runs,
    }
}

fn encode(ops: &[u8]) -> (String, Vec<CigarOp>) {
    let mut out = String::new();
    let mut runs = Vec::new();
    let mut k = 0;
    while k < ops.len() {
        let c = ops[k];
        let start = k;
        while k < ops.len() && ops[k] == c {
            k += 1;
        }
        out.push_str(&(k - start).to_string());
        out.push(c as char);
        runs.push((k - start, c));
    }
    (out, runs)
}

#[cfg(test)]
mod tests {

    /// Row order vs wavefront on realistic consensus lengths, so the choice is
    /// made on a measurement rather than on which one looks faster.
    ///
    /// Not a `#[test]` assertion --- it prints and is read. Run with
    /// `cargo test --release --lib bench_fill_orders -- --nocapture --ignored`.
    #[test]
    #[ignore]
    fn bench_fill_orders() {
        let mut seed = 0x9E3779B97F4A7C15u64;
        let mut next = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        // 1364 bp is the measured mean consensus length on droso_deep.
        let len = 1364usize;
        let pairs: Vec<(Vec<u8>, Vec<u8>)> = (0..12)
            .map(|_| {
                let a: Vec<u8> = (0..len).map(|_| b"ACGT"[(next() % 4) as usize]).collect();
                // A near-copy, which is what merge candidates actually look like.
                let mut b = a.clone();
                for _ in 0..len / 50 {
                    let p = (next() as usize) % len;
                    b[p] = b"ACGT"[(next() % 4) as usize];
                }
                (a, b)
            })
            .collect();

        // The full call, so the fill's share of it is visible: vectorising the
        // fill can only ever win back the fill's share.
        let t_full = std::time::Instant::now();
        let mut sink2 = 0i64;
        for (a, b) in &pairs {
            sink2 += semiglobal_with(a, b, Scoring::MERGE, TieBreak::PARASAIL).score as i64;
        }
        eprintln!(
            "  FULL semiglobal_with (fill + traceback + strings): {:.3}s  (sink {sink2})",
            t_full.elapsed().as_secs_f64()
        );

        for (name, f) in [("row      ", 0u8), ("wavefront", 1u8)] {
            let t0 = std::time::Instant::now();
            let mut sink = 0i64;
            for (a, b) in &pairs {
                let t = if f == 0 {
                    forward(a, b, Scoring::MERGE, TieBreak::PARASAIL, None, false)
                } else {
                    forward_wavefront(a, b, Scoring::MERGE, TieBreak::PARASAIL)
                };
                sink += t.last_row[t.m] as i64;
            }
            let el = t0.elapsed().as_secs_f64();
            let cells = pairs.len() as f64 * (len * len) as f64;
            eprintln!(
                "  {name}: {:.2}s for {} pairs of {len}bp -> {:.0}M cells/s (sink {sink})",
                el,
                pairs.len(),
                cells / el / 1e6
            );
        }
    }

    /// Global alignment must *charge* end gaps where semi-global does not.
    /// This is the property the bubble-popping change depends on, so it is
    /// pinned rather than assumed.
    #[test]
    fn global_charges_end_gaps_and_semiglobal_does_not() {
        let core = b"ACGTTGCAAGGCTTAACCGGATTCAGGTACGATCGATCGGCTAACGT".to_vec();
        let mut with_tail = core.clone();
        with_tail.extend_from_slice(b"TTTTTTTTTTTTTTTTTTTT"); // 20-base overhang

        let sg_same = semiglobal(&core, &core, Scoring::BUBBLE).score;
        let sg_tail = semiglobal(&with_tail, &core, Scoring::BUBBLE).score;
        assert_eq!(
            sg_tail, sg_same,
            "semi-global: a trailing overhang is free, so the score is unchanged"
        );

        let gl_same = global(&core, &core, Scoring::BUBBLE).score;
        let gl_tail = global(&with_tail, &core, Scoring::BUBBLE).score;
        assert_eq!(gl_same, sg_same, "with no overhang the two agree");
        assert!(
            gl_tail < gl_same,
            "global: the overhang must cost something ({gl_tail} vs {gl_same})"
        );
        // open 12 + 19 extensions at 1 = 31.
        assert_eq!(gl_same - gl_tail, 31, "one gap open plus 19 extensions");
    }

    #[test]
    fn global_spans_both_sequences_completely() {
        let a = b"ACGTACGTAAGGCCTTACGT".to_vec();
        let b = b"ACGTACGTTAGGCCTTACGT".to_vec();
        let aln = global(&a, &b, Scoring::BUBBLE);
        let (mut qi, mut rj) = (0usize, 0usize);
        for &(l, op) in &aln.ops {
            match op {
                b'=' | b'X' | b'M' => {
                    qi += l;
                    rj += l;
                }
                b'I' => qi += l,
                b'D' => rj += l,
                _ => {}
            }
        }
        assert_eq!((qi, rj), (a.len(), b.len()), "every base is accounted for");
    }

    /// A tie-break deliberately different from parasail's, so the comparison
    /// covers the `prefer_open` and ordering branches rather than one path.
    fn alt_tiebreak() -> TieBreak {
        let mut tb = TieBreak::PARASAIL;
        tb.prefer_open = !tb.prefer_open;
        tb.h_order = [Step::Insert, Step::Delete, Step::Diagonal];
        tb
    }

    /// The anti-diagonal fill must reproduce the row fill exactly --- every
    /// traceback byte, every last-row and last-column score. This is what makes
    /// the wavefront a safe base for the SIMD kernel: if this holds, the only
    /// thing SIMD has to match is a function already known to equal `forward`.
    #[test]
    fn wavefront_matches_the_row_order_exactly() {
        let mut seed = 0x2545F4914F6CDD1Du64;
        let mut next = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        let mut checked = 0usize;
        for case in 0..400 {
            // Include the degenerate shapes: empty, length 1, and very uneven.
            let n = (next() % 40) as usize + case % 3;
            let m = (next() % 40) as usize + (case / 3) % 3;
            let mk = |len: usize, f: &mut dyn FnMut() -> u64| -> Vec<u8> {
                (0..len).map(|_| b"ACGT"[(f() % 4) as usize]).collect()
            };
            let s1 = mk(n, &mut next);
            let s2 = mk(m, &mut next);
            for tb in [TieBreak::PARASAIL, alt_tiebreak()] {
                let a = forward(&s1, &s2, Scoring::MERGE, tb, None, false);
                let b = forward_wavefront(&s1, &s2, Scoring::MERGE, tb);
                assert_eq!(a.packed, b.packed, "traceback bytes differ, n={n} m={m}");
                assert_eq!(a.last_col, b.last_col, "last column differs, n={n} m={m}");
                assert_eq!(a.last_row, b.last_row, "last row differs, n={n} m={m}");
                checked += 1;
            }
        }
        assert!(checked >= 800, "only {checked} comparisons ran");
    }
    use super::*;

    /// The scoring these semantics tests were originally verified against.
    ///
    /// `match=4, mismatch=-8, open=12, ext=1` --- isONcorrect's, not isONform's.
    /// Kept deliberately: the expected scores below were read off the parasail
    /// *library* at this scoring, so they are measurements, and re-deriving them
    /// arithmetically under `Scoring::BUBBLE` would turn evidence into
    /// assumption. What they establish --- free end gaps on both sequences, a
    /// one-base gap costing `open`, longer gaps costing `open + (L-1) * ext` ---
    /// is a property of parasail's semi-global mode and holds at any scoring.
    ///
    /// isONform's own scorings are [`Scoring::BUBBLE`] and [`Scoring::MERGE`];
    /// the tie-break and CIGAR tests below use those.
    const VERIFIED: Scoring = Scoring {
        match_score: 4,
        mismatch: -8,
        open: 12,
        ext: 1,
    };

    fn sg(a: &str, b: &str) -> Alignment {
        semiglobal(a.as_bytes(), b.as_bytes(), VERIFIED)
    }

    /// `parasail_matrix_create("ACGT", m, n)` scores anything outside `ACGT` as
    /// **zero against everything, itself included**. Every number below was read
    /// off the parasail library, not derived:
    ///
    /// ```text
    /// m = parasail.matrix_create("ACGT", 2, -8)
    /// parasail.sg_trace_scan_16("A", "A", 12, 1, m).score  ->  2
    /// parasail.sg_trace_scan_16("A", "X", 12, 1, m).score  ->  0
    /// parasail.sg_trace_scan_16("X", "X", 12, 1, m).score  ->  0
    /// parasail.sg_trace_scan_16("N", "N", 12, 1, m).score  ->  0
    /// parasail.sg_trace_scan_16("a", "A", 12, 1, m).score  ->  2   (case-insensitive)
    /// ```
    #[test]
    fn a_character_outside_acgt_scores_zero_even_against_itself() {
        let b = Scoring::BUBBLE;
        assert_eq!(subst(b, b'A', b'A'), 2, "a base against itself matches");
        assert_eq!(subst(b, b'A', b'C'), -8, "two bases mismatch");
        assert_eq!(subst(b, b'A', b'X'), 0, "base against non-base is zero");
        assert_eq!(subst(b, b'X', b'A'), 0, "and symmetrically");
        assert_eq!(
            subst(b, b'X', b'X'),
            0,
            "X against X is ZERO, not a match --- this is the whole point"
        );
        assert_eq!(subst(b, b'N', b'N'), 0);
        assert_eq!(subst(b, b'X', b'N'), 0);
        // parasail's alphabet lookup is case-insensitive.
        assert_eq!(subst(b, b'a', b'A'), 2);
        assert_eq!(subst(b, b'g', b'g'), 2);
        assert_eq!(subst(b, b'a', b'c'), -8);
    }

    /// The case that found finding 24's scoring bug, and later its end-cell bug.
    /// `generate_consensus_path` returns `"X" * max_len` when every path span in a
    /// bubble is shorter than `k`, so a 17-X against an 18-X placeholder is a real
    /// comparison the simplification stage makes.
    ///
    /// parasail, verbatim:
    ///
    /// ```text
    /// m = parasail.matrix_create("ACGT", 2, -8)
    /// r = parasail.sg_trace_scan_16("X"*17, "X"*18, 12, 1, m)
    /// r.score                      ->  0
    /// str(r.cigar.decode, "utf-8") ->  "16I1=17D"
    /// ```
    ///
    /// Both halves now match. The score needed `subst` (every X pair scores 0,
    /// not `match_score`); the CIGAR needed the end-cell scan to start at 1, so
    /// that "align nothing" is not a candidate. Before the second fix this
    /// asserted the port's own `17I18D` as a known divergence — a residual that
    /// looked structural and was not.
    #[test]
    fn all_x_placeholders_match_parasail_exactly() {
        let a = "X".repeat(17);
        let b = "X".repeat(18);
        let got = semiglobal(a.as_bytes(), b.as_bytes(), Scoring::BUBBLE);
        assert_eq!(got.score, 0, "no pairing of X's can beat any other");
        assert_eq!(got.cigar, "16I1=17D", "parasail's own CIGAR for this input");
    }

    /// parasail's `sg` insists on at least one diagonal step: it will not take the
    /// free all-end-gap path even when every real alignment scores below zero.
    /// Read off the library at `match=2, mismatch=-2, open=12, ext=1`:
    ///
    /// ```text
    /// sg("A",    "C")    -> -2  "1X"
    /// sg("AA",   "CC")   -> -2  "1I1X1D"
    /// sg("AAAA", "CCCC") -> -2  "3I1X3D"
    /// sg("A",    "CCCC") -> -2  "1X3D"
    /// sg("AC",   "CA")   ->  2  "1I1=1D"      one pair does match
    /// ```
    ///
    /// A port that admits end cell `(n, 0)` or `(0, m)` returns 0 on every one of
    /// the −2 rows, because walking the free border costs nothing. That was 12 of
    /// 54 884 recorded real calls.
    #[test]
    fn sg_requires_at_least_one_diagonal_rather_than_taking_a_free_zero() {
        let sc = Scoring {
            match_score: 2,
            mismatch: -2,
            open: 12,
            ext: 1,
        };
        // Score and CIGAR both, read off the library. `("AAAA", "C")` is the one
        // that pins the corner rule: with `m == 1` the last-row scan is empty, so
        // the choice falls to the last column and lands on `(1, 1)` --- which is
        // why the CIGAR is `1X3I` and not `3I1X`.
        let cases = [
            ("A", "C", -2, "1X"),
            ("AA", "CC", -2, "1I1X1D"),
            ("AAA", "CCC", -2, "2I1X2D"),
            ("A", "CCCC", -2, "1X3D"),
            ("AAAA", "C", -2, "1X3I"),
            ("AAAA", "CCCC", -2, "3I1X3D"),
            ("AC", "CA", 2, "1I1=1D"),
            ("ACGT", "TGCA", 2, "3I1=3D"),
        ];
        for (a, b, score, cigar) in cases {
            let got = semiglobal(a.as_bytes(), b.as_bytes(), sc);
            assert_eq!(got.score, score, "score for sg({a:?}, {b:?})");
            assert_eq!(got.cigar, cigar, "cigar for sg({a:?}, {b:?})");
        }
    }

    #[test]
    fn the_mismatch_penalty_changes_where_a_gap_beats_a_mismatch() {
        // Same input, the two scorings, different alignments --- which is why they
        // are separate constants rather than one shared default.
        let a = b"ACGTACGT";
        let b = b"ACGAACGT";
        let strict = semiglobal(a, b, Scoring::BUBBLE);
        let loose = semiglobal(a, b, Scoring::MERGE);
        assert_ne!(
            strict.score, loose.score,
            "a mismatch is priced differently by the two"
        );
    }

    /// Ground truth read off the parasail library itself, not assumed.
    #[test]
    fn end_gaps_are_free_on_both_sequences() {
        let a = sg("ACGTACGT", "TTTACGTACGTTTT");
        assert_eq!(a.cigar, "3D8=3D");
        assert_eq!(a.score, 32);

        let b = sg("TTTACGTACGTTTT", "ACGTACGT");
        assert_eq!(b.cigar, "3I8=3I");
        assert_eq!(b.score, 32);
    }

    #[test]
    fn identical_sequences_score_every_match() {
        let a = sg("ACGTACGT", "ACGTACGT");
        assert_eq!(a.cigar, "8=");
        assert_eq!(a.score, 32);
    }

    #[test]
    fn a_single_base_gap_costs_the_opening_penalty() {
        // Both sides long enough that end gaps cannot win.
        let s1 = "ACGTACGTACGTTGCAACGT";
        let s2 = "ACGTACGTACGTGCAACGT"; // one T dropped
        let a = sg(s1, s2);
        assert_eq!(a.score, 19 * 4 - 12);
        assert_eq!(a.cigar.matches('I').count(), 1);
    }

    #[test]
    fn a_longer_gap_costs_open_plus_extensions() {
        let s1 = "ACGTACGTACGTTTTTGCAACGTACGT";
        let s2 = "ACGTACGTACGTGCAACGTACGT"; // four T dropped
        let a = sg(s1, s2);
        assert_eq!(a.score, 23 * 4 - (12 + 3));
    }

    #[test]
    fn the_cigar_spans_both_sequences() {
        for (s1, s2) in [
            ("ACGTACGT", "TTTACGTACGTTTT"),
            ("ACGTACGTACGTTGCAACGT", "ACGTACGTACGTGCAACGT"),
            ("ACGT", "ACGT"),
            ("", "ACGT"),
            ("ACGT", ""),
        ] {
            let a = semiglobal(s1.as_bytes(), s2.as_bytes(), Scoring::BUBBLE);
            let (mut q, mut r) = (0usize, 0usize);
            for (len, op) in &a.ops {
                match op {
                    b'=' | b'X' => {
                        q += len;
                        r += len;
                    }
                    b'I' => q += len,
                    b'D' => r += len,
                    _ => panic!("unexpected op"),
                }
            }
            assert_eq!((q, r), (s1.len(), s2.len()), "on {s1:?} vs {s2:?}");
        }
    }

    #[test]
    fn an_empty_sequence_is_all_free_gaps() {
        let a = sg("", "ACGT");
        assert_eq!((a.cigar.as_str(), a.score), ("4D", 0));
        let b = sg("ACGT", "");
        assert_eq!((b.cigar.as_str(), b.score), ("4I", 0));
        let c = sg("", "");
        assert_eq!((c.cigar.as_str(), c.score), ("", 0));
    }
}

/// Differential oracle against the parasail library.
///
/// ```bash
/// bench/dump_reference.py --fastq cluster.fastq --outdir /tmp/d
/// PARASAIL_CASES=/tmp/d/parasail_cases.tsv \
///   cargo test --manifest-path rust/Cargo.toml parasail::oracle -- --nocapture
/// ```
///
/// `sweep` is the tool that *chose* [`TieBreak::PARASAIL`]: it runs every
/// combination against the recorded cases and prints the score. Re-run it if a
/// mismatch ever appears — a new corpus reaching a tie the current ordering gets
/// wrong would show up as a different winner.
#[cfg(test)]
mod oracle {
    use super::*;

    struct Case {
        s1: Vec<u8>,
        s2: Vec<u8>,
        cigar: String,
        score: i32,
        scoring: Scoring,
    }

    fn load() -> Option<Vec<Case>> {
        let path = std::env::var("PARASAIL_CASES").ok()?;
        let data = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("PARASAIL_CASES={path} unreadable: {e}"));
        let mut cases = Vec::new();
        for line in data
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let f: Vec<&str> = line.split('\t').collect();
            assert!(f.len() >= 8, "short row");
            cases.push(Case {
                s1: f[0].as_bytes().to_vec(),
                s2: f[1].as_bytes().to_vec(),
                cigar: f[2].to_string(),
                score: f[3].parse().unwrap(),
                scoring: Scoring {
                    match_score: f[4].parse().unwrap(),
                    mismatch: f[5].parse().unwrap(),
                    open: f[6].parse().unwrap(),
                    ext: f[7].parse().unwrap(),
                },
            });
        }
        Some(cases)
    }

    #[test]
    fn matches_parasail_exactly() {
        let Some(cases) = load() else { return };
        assert!(!cases.is_empty(), "no cases");

        let (mut bad_score, mut bad_cigar, mut shown) = (0usize, 0usize, 0usize);
        for c in &cases {
            let got = semiglobal(&c.s1, &c.s2, c.scoring);
            if got.score != c.score {
                if bad_score < 5 {
                    eprintln!(
                        "SCOREMISMATCH len {}x{} scoring m={} n={} o={} e={} parasail={} rust={}\n  s1={}\n  s2={}",
                        c.s1.len(),
                        c.s2.len(),
                        c.scoring.match_score,
                        c.scoring.mismatch,
                        c.scoring.open,
                        c.scoring.ext,
                        c.score,
                        got.score,
                        String::from_utf8_lossy(&c.s1),
                        String::from_utf8_lossy(&c.s2),
                    );
                }
                bad_score += 1;
            }
            if got.cigar != c.cigar {
                if shown < 3 {
                    shown += 1;
                    eprintln!(
                        "cigar mismatch on {}x{}:\n  parasail: {}\n  rust    : {}",
                        c.s1.len(),
                        c.s2.len(),
                        c.cigar,
                        got.cigar
                    );
                }
                bad_cigar += 1;
            }
        }
        eprintln!(
            "parasail oracle: {} cases, {bad_score} score mismatches, {bad_cigar} cigar mismatches",
            cases.len()
        );
        // The **score** is exact and gated as such. It was not always: admitting
        // end cell `(n, 0)` / `(0, m)` returned 0 where parasail returns a
        // negative optimum, on 12 of 54 884 recorded real calls. See `end_cell`.
        assert_eq!(bad_score, 0, "scores differed");

        // The **CIGAR** is exact too, and gated as such. It was not: parasail
        // resolves a tie among equally-scoring end cells by excluding the corner
        // from the last-row and last-column scans and considering it last, which
        // `end_cell` now models. Before that, 112 of these 54 884 differed --- all
        // of them equally *optimal* paths, which is why the score gate above could
        // not see them.
        assert_eq!(
            bad_cigar, 0,
            "CIGARs differed --- if these are equally-optimal alternatives rather \
             than worse paths, `rescore` below is the thing to reach for"
        );

        // Kept as a diagnostic for the next time a CIGAR does differ: it separates
        // "picked a different member of a tie" from "picked a worse path", which
        // the score gate alone cannot do, since a tie has the same score by
        // definition.
        for c in &cases {
            let got = semiglobal(&c.s1, &c.s2, c.scoring);
            debug_assert_eq!(
                rescore(&got.ops, &c.s1, &c.s2, c.scoring),
                c.score,
                "the port's own CIGAR does not rescore to its reported score"
            );
        }
    }

    /// Score a CIGAR the way parasail's semi-global mode does: leading and
    /// trailing gaps are free, interior gaps cost `open + (L-1) * ext`, and
    /// diagonals use the same substitution rule as the aligner.
    ///
    /// Used to check that a CIGAR the port reports is genuinely optimal rather
    /// than merely different, which is the one thing a differing CIGAR must still
    /// satisfy.
    fn rescore(ops: &[CigarOp], s1: &[u8], s2: &[u8], sc: Scoring) -> i32 {
        let mut score = 0i32;
        let (mut i, mut j) = (0usize, 0usize);
        for (idx, &(len, op)) in ops.iter().enumerate() {
            let terminal = idx == 0 || idx == ops.len() - 1;
            match op {
                b'=' | b'X' | b'M' => {
                    for _ in 0..len {
                        score += subst(sc, s1[i], s2[j]);
                        i += 1;
                        j += 1;
                    }
                }
                // `I` consumes s1, `D` consumes s2 --- see the module docs.
                b'I' | b'D' => {
                    if !terminal {
                        score -= sc.open + (len as i32 - 1) * sc.ext;
                    }
                    if op == b'I' {
                        i += len;
                    } else {
                        j += len;
                    }
                }
                other => panic!("unexpected CIGAR op {}", other as char),
            }
        }
        assert_eq!(i, s1.len(), "CIGAR does not consume all of s1");
        assert_eq!(j, s2.len(), "CIGAR does not consume all of s2");
        score
    }

    /// Writes only the cases the current settings get wrong, so the sweep can run
    /// on the discriminating subset instead of all 54 884.
    #[test]
    fn dump_hard_cases() {
        let Some(cases) = load() else { return };
        let Ok(out) = std::env::var("PARASAIL_HARD_OUT") else {
            return;
        };
        let mut w = String::new();
        for c in &cases {
            let got = semiglobal(&c.s1, &c.s2, c.scoring);
            if got.score != c.score || got.cigar != c.cigar {
                w.push_str(&format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    String::from_utf8_lossy(&c.s1),
                    String::from_utf8_lossy(&c.s2),
                    c.cigar,
                    c.score,
                    c.scoring.match_score,
                    c.scoring.mismatch,
                    c.scoring.open,
                    c.scoring.ext,
                    got.cigar,
                    got.score
                ));
            }
        }
        std::fs::write(&out, w).unwrap();
        eprintln!("wrote hard cases to {out}");
    }

    /// Not an assertion --- a measurement. Prints how each tie-break ordering
    /// scores so the right one can be picked from evidence.
    #[test]
    fn sweep() {
        let Some(mut cases) = load() else { return };
        let Ok(setting) = std::env::var("PARASAIL_SWEEP") else {
            return;
        };
        // 96 configurations times full-length reads is hours of quadratic DP, so
        // the variable doubles as a case cap: PARASAIL_SWEEP=200 sweeps the
        // first 200. Ties are what matter here, not volume.
        if let Ok(cap) = setting.parse::<usize>() {
            if cap > 1 {
                cases.truncate(cap);
            }
        }
        let mut results: Vec<(usize, TieBreak)> = TieBreak::all()
            .into_iter()
            .map(|tb| {
                let hits = cases
                    .iter()
                    .filter(|c| semiglobal_with(&c.s1, &c.s2, c.scoring, tb).cigar == c.cigar)
                    .count();
                (hits, tb)
            })
            .collect();
        results.sort_by_key(|(hits, _)| std::cmp::Reverse(*hits));
        eprintln!("tie-break sweep over {} cases:", cases.len());
        for (hits, tb) in results.iter().take(8) {
            eprintln!(
                "  {hits:>6} / {}  h_order={:?} open={} col_first={} last_max={}",
                cases.len(),
                tb.h_order,
                tb.prefer_open,
                tb.column_first,
                tb.last_max
            );
        }
    }
}

/// Does banding change the answer? A measurement, not a feature.
///
/// The guard is ~65% of runtime and its DP is `O(n*m)` over a read against its
/// own correction. Those two sequences differ very little, so a band around the
/// diagonal should lose nothing --- but `fix_correction` walks the *path*, so an
/// equally-optimal-but-different path changes the output. This computes the
/// banded answer and compares CIGARs against what parasail recorded.
///
/// ```bash
/// PARASAIL_CASES=/tmp/d/parasail_cases.tsv \
///   cargo test --release --manifest-path rust/Cargo.toml parasail::banding -- --nocapture
/// ```
#[cfg(test)]
mod banding {
    use super::*;

    /// Same recurrence, restricted to `|i - j - drift| <= half`, where `drift`
    /// centres the band on the diagonal the length difference implies.
    fn banded(s1: &[u8], s2: &[u8], sc: Scoring, tb: TieBreak, half: usize) -> Option<String> {
        let (n, m) = (s1.len(), s2.len());
        if n == 0 || m == 0 {
            return Some(semiglobal_with(s1, s2, sc, tb).cigar);
        }
        let lo_of = |i: usize| -> usize { i.saturating_sub(half).max(1) };
        let hi_of = |i: usize| -> usize { (i + half).min(m) };

        let mut h_prev = vec![0i32; m + 1];
        let mut e_prev = vec![0i32; m + 1];
        let mut f_prev = vec![NEG; m + 1];
        f_prev[0] = 0;
        let mut h_cur = vec![NEG; m + 1];
        let mut e_cur = vec![NEG; m + 1];
        let mut f_cur = vec![NEG; m + 1];

        let mut packed = vec![0u8; n * m];
        let mut last_col = vec![NEG; n + 1];
        last_col[0] = 0;
        let mut touched_edge = false;

        for i in 1..=n {
            let (lo, hi) = (lo_of(i), hi_of(i));
            // Cells outside the band are unreachable, not zero.
            for v in h_cur.iter_mut() {
                *v = NEG;
            }
            for v in e_cur.iter_mut() {
                *v = NEG;
            }
            for v in f_cur.iter_mut() {
                *v = NEG;
            }
            if lo == 1 {
                h_cur[0] = 0;
                f_cur[0] = 0;
            }
            let c1 = s1[i - 1];
            for j in lo..=hi {
                let open_d = h_cur[j - 1] - sc.open;
                let ext_d = e_cur[j - 1] - sc.ext;
                e_cur[j] = open_d.max(ext_d);
                let d_leaves = if tb.prefer_open {
                    open_d >= ext_d
                } else {
                    open_d > ext_d
                };
                let open_i = h_prev[j] - sc.open;
                let ext_i = f_prev[j] - sc.ext;
                f_cur[j] = open_i.max(ext_i);
                let i_leaves = if tb.prefer_open {
                    open_i >= ext_i
                } else {
                    open_i > ext_i
                };
                let diag = h_prev[j - 1] + subst(sc, c1, s2[j - 1]);
                let best = diag.max(e_cur[j]).max(f_cur[j]);
                h_cur[j] = best;
                let mut choice = 0u8;
                for step in tb.h_order {
                    let (hit, code) = match step {
                        Step::Diagonal => (diag == best, 0u8),
                        Step::Delete => (e_cur[j] == best, 1),
                        Step::Insert => (f_cur[j] == best, 2),
                    };
                    if hit {
                        choice = code;
                        break;
                    }
                }
                packed[(i - 1) * m + (j - 1)] =
                    choice | (u8::from(d_leaves) << 2) | (u8::from(i_leaves) << 3);
            }
            last_col[i] = if hi == m { h_cur[m] } else { NEG };
            std::mem::swap(&mut h_prev, &mut h_cur);
            std::mem::swap(&mut e_prev, &mut e_cur);
            std::mem::swap(&mut f_prev, &mut f_cur);
        }

        let mut best = (NEG, n, m);
        for (j, &v) in h_prev.iter().enumerate() {
            if v > best.0 {
                best = (v, n, j);
            }
        }
        for (i, &v) in last_col.iter().enumerate() {
            if v > best.0 {
                best = (v, i, m);
            }
        }
        let (_score, mut i, mut j) = best;

        let mut ops: Vec<u8> = Vec::with_capacity(n + m);
        ops.extend(std::iter::repeat_n(b'D', m - j));
        ops.extend(std::iter::repeat_n(b'I', n - i));
        let mut chain: Option<Step> = None;
        while i > 0 && j > 0 {
            // A *real* band boundary, not the sequence boundary: `lo`/`hi` are
            // clamped to [1, m], and hitting those is not clipping.
            let (lo, hi) = (lo_of(i), hi_of(i));
            if (j == lo && lo > 1) || (j == hi && hi < m) {
                touched_edge = true;
            }
            let bits = packed[(i - 1) * m + (j - 1)];
            match chain {
                Some(Step::Delete) => {
                    ops.push(b'D');
                    j -= 1;
                    chain = if bits & 0b100 != 0 {
                        None
                    } else {
                        Some(Step::Delete)
                    };
                }
                Some(Step::Insert) => {
                    ops.push(b'I');
                    i -= 1;
                    chain = if bits & 0b1000 != 0 {
                        None
                    } else {
                        Some(Step::Insert)
                    };
                }
                Some(Step::Diagonal) => unreachable!(),
                None => match bits & 0b11 {
                    0 => {
                        ops.push(if s1[i - 1] == s2[j - 1] { b'=' } else { b'X' });
                        i -= 1;
                        j -= 1;
                    }
                    1 => chain = Some(Step::Delete),
                    _ => chain = Some(Step::Insert),
                },
            }
        }
        ops.extend(std::iter::repeat_n(b'D', j));
        ops.extend(std::iter::repeat_n(b'I', i));
        ops.reverse();
        if touched_edge {
            return None; // the band clipped the path; a wider one is needed
        }
        Some(encode(&ops).0)
    }

    /// Would "same score at W and at 2W" have caught the cases a too-narrow
    /// band gets wrong? If so it is a sound-in-practice acceptance test: a
    /// better path within 2W is found by the wider run and the scores diverge.
    #[test]
    fn does_widening_detect_a_band_that_is_too_narrow() {
        let Ok(path) = std::env::var("PARASAIL_CASES") else {
            return;
        };
        let data = std::fs::read_to_string(&path).expect("cases");
        let (mut wrong_caught, mut wrong_missed, mut ok) = (0usize, 0usize, 0usize);
        for line in data
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 8 {
                continue;
            }
            let (s1, s2) = (f[0].as_bytes(), f[1].as_bytes());
            let sc = Scoring {
                match_score: f[4].parse().unwrap(),
                mismatch: f[5].parse().unwrap(),
                open: f[6].parse().unwrap(),
                ext: f[7].parse().unwrap(),
            };
            // Deliberately too narrow, so wrong answers actually occur.
            let narrow = banded_score(s1, s2, sc, TieBreak::PARASAIL, 32);
            let wider = banded_score(s1, s2, sc, TieBreak::PARASAIL, 64);
            let truth: i32 = f[3].parse().unwrap();
            match (narrow, wider) {
                (Some(a), Some(b)) => {
                    if a == truth {
                        ok += 1;
                    } else if a != b {
                        wrong_caught += 1;
                    } else {
                        wrong_missed += 1;
                    }
                }
                _ => ok += 1,
            }
        }
        eprintln!(
            "narrow band correct {ok}, wrong-and-caught-by-widening {wrong_caught}, \
             wrong-and-missed {wrong_missed}"
        );
    }

    /// Banded score only, without the traceback or the edge check.
    fn banded_score(s1: &[u8], s2: &[u8], sc: Scoring, tb: TieBreak, half: usize) -> Option<i32> {
        let (n, m) = (s1.len(), s2.len());
        if n == 0 || m == 0 {
            return None;
        }
        let mut h_prev = vec![0i32; m + 1];
        let mut e_prev = vec![0i32; m + 1];
        let mut f_prev = vec![NEG; m + 1];
        f_prev[0] = 0;
        let mut h_cur = vec![NEG; m + 1];
        let mut e_cur = vec![NEG; m + 1];
        let mut f_cur = vec![NEG; m + 1];
        let mut best = NEG;
        for i in 1..=n {
            let lo = i.saturating_sub(half).max(1);
            let hi = (i + half).min(m);
            h_cur.fill(NEG);
            e_cur.fill(NEG);
            f_cur.fill(NEG);
            if lo == 1 {
                h_cur[0] = 0;
                f_cur[0] = 0;
            }
            let c1 = s1[i - 1];
            for j in lo..=hi {
                e_cur[j] = (h_cur[j - 1] - sc.open).max(e_cur[j - 1] - sc.ext);
                f_cur[j] = (h_prev[j] - sc.open).max(f_prev[j] - sc.ext);
                let diag = h_prev[j - 1] + subst(sc, c1, s2[j - 1]);
                h_cur[j] = diag.max(e_cur[j]).max(f_cur[j]);
            }
            if hi == m {
                best = best.max(h_cur[m]);
            }
            std::mem::swap(&mut h_prev, &mut h_cur);
            std::mem::swap(&mut e_prev, &mut e_cur);
            std::mem::swap(&mut f_prev, &mut f_cur);
            let _ = tb;
        }
        for &v in h_prev.iter() {
            best = best.max(v);
        }
        Some(best)
    }

    #[test]
    fn how_often_does_banding_change_the_cigar() {
        let Ok(path) = std::env::var("PARASAIL_CASES") else {
            return;
        };
        let data = std::fs::read_to_string(&path).expect("cases");
        let mut cases = Vec::new();
        for line in data
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 8 {
                continue;
            }
            cases.push((
                f[0].as_bytes().to_vec(),
                f[1].as_bytes().to_vec(),
                f[2].to_string(),
                Scoring {
                    match_score: f[4].parse().unwrap(),
                    mismatch: f[5].parse().unwrap(),
                    open: f[6].parse().unwrap(),
                    ext: f[7].parse().unwrap(),
                },
            ));
        }
        eprintln!("{} recorded alignments", cases.len());
        for half in [32usize, 64, 128, 256, 512] {
            let (mut same, mut differ, mut clipped) = (0usize, 0usize, 0usize);
            for (s1, s2, want, sc) in &cases {
                match banded(s1, s2, *sc, TieBreak::PARASAIL, half) {
                    None => clipped += 1,
                    Some(got) if &got == want => same += 1,
                    Some(_) => differ += 1,
                }
            }
            eprintln!(
                "  half-band {half:>4}: identical {same:>5}  DIFFERENT {differ:>5}  \
                 clipped(fallback) {clipped:>5}"
            );
        }
    }
}

/// Does the guard actually need *semi*-global alignment?
///
/// `correct_read` aligns `seq` against `corr`, and `corr` is built by splicing
/// corrected spans into `seq` — so the two share an **exact prefix and suffix**.
/// Matching those scores `+match` per base while a free end gap scores 0, so
/// free ends should never win. If exact affine *global* alignment agrees with
/// parasail's `sg` on real data, the semi-globalness is an artifact of reusing
/// `parasail_alignment`'s default, and any global aligner is a candidate.
///
/// ```bash
/// PARASAIL_CASES=/tmp/d/parasail_cases.tsv \
///   cargo test --release --manifest-path rust/Cargo.toml parasail::global_vs_sg -- --nocapture
/// ```
#[cfg(test)]
mod global_vs_sg {
    use super::*;

    /// Exact affine global (Needleman-Wunsch-Gotoh), same tie-break as
    /// [`TieBreak::PARASAIL`]. Self-contained so the production path is
    /// untouched by this measurement.
    fn global(s1: &[u8], s2: &[u8], sc: Scoring, tb: TieBreak) -> String {
        let (n, m) = (s1.len(), s2.len());
        let at = |i: usize, j: usize| i * (m + 1) + j;
        let mut h = vec![NEG; (n + 1) * (m + 1)];
        let mut e = vec![NEG; (n + 1) * (m + 1)];
        let mut f = vec![NEG; (n + 1) * (m + 1)];
        let gap = |l: usize| -> i32 {
            if l == 0 {
                0
            } else {
                -(sc.open + (l as i32 - 1) * sc.ext)
            }
        };
        h[at(0, 0)] = 0;
        for i in 1..=n {
            h[at(i, 0)] = gap(i);
            f[at(i, 0)] = gap(i);
        }
        for j in 1..=m {
            h[at(0, j)] = gap(j);
            e[at(0, j)] = gap(j);
        }
        for i in 1..=n {
            for j in 1..=m {
                e[at(i, j)] = (h[at(i, j - 1)] - sc.open).max(e[at(i, j - 1)] - sc.ext);
                f[at(i, j)] = (h[at(i - 1, j)] - sc.open).max(f[at(i - 1, j)] - sc.ext);
                let diag = h[at(i - 1, j - 1)] + subst(sc, s1[i - 1], s2[j - 1]);
                h[at(i, j)] = diag.max(e[at(i, j)]).max(f[at(i, j)]);
            }
        }

        let (mut i, mut j) = (n, m);
        let mut ops: Vec<u8> = Vec::new();
        let mut chain: Option<Step> = None;
        while i > 0 && j > 0 {
            match chain {
                Some(Step::Delete) => {
                    let open_d = h[at(i, j - 1)] - sc.open;
                    let ext_d = e[at(i, j - 1)] - sc.ext;
                    let leaves = if tb.prefer_open {
                        open_d >= ext_d
                    } else {
                        open_d > ext_d
                    };
                    ops.push(b'D');
                    j -= 1;
                    chain = if leaves { None } else { Some(Step::Delete) };
                }
                Some(Step::Insert) => {
                    let open_i = h[at(i - 1, j)] - sc.open;
                    let ext_i = f[at(i - 1, j)] - sc.ext;
                    let leaves = if tb.prefer_open {
                        open_i >= ext_i
                    } else {
                        open_i > ext_i
                    };
                    ops.push(b'I');
                    i -= 1;
                    chain = if leaves { None } else { Some(Step::Insert) };
                }
                Some(Step::Diagonal) => unreachable!(),
                None => {
                    let cur = h[at(i, j)];
                    let same = s1[i - 1] == s2[j - 1];
                    let diag = h[at(i - 1, j - 1)] + subst(sc, s1[i - 1], s2[j - 1]);
                    let mut taken = Step::Diagonal;
                    for step in tb.h_order {
                        let hit = match step {
                            Step::Diagonal => diag == cur,
                            Step::Delete => e[at(i, j)] == cur,
                            Step::Insert => f[at(i, j)] == cur,
                        };
                        if hit {
                            taken = step;
                            break;
                        }
                    }
                    match taken {
                        Step::Diagonal => {
                            ops.push(if same { b'=' } else { b'X' });
                            i -= 1;
                            j -= 1;
                        }
                        step => chain = Some(step),
                    }
                }
            }
        }
        ops.extend(std::iter::repeat_n(b'D', j));
        ops.extend(std::iter::repeat_n(b'I', i));
        ops.reverse();
        encode(&ops).0
    }

    #[test]
    fn global_agrees_with_parasail_semiglobal() {
        let Ok(path) = std::env::var("PARASAIL_CASES") else {
            return;
        };
        let data = std::fs::read_to_string(&path).expect("cases");
        let (mut n, mut same, mut differ, mut shown) = (0usize, 0usize, 0usize, 0usize);
        for line in data
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 8 {
                continue;
            }
            let sc = Scoring {
                match_score: f[4].parse().unwrap(),
                mismatch: f[5].parse().unwrap(),
                open: f[6].parse().unwrap(),
                ext: f[7].parse().unwrap(),
            };
            let got = global(f[0].as_bytes(), f[1].as_bytes(), sc, TieBreak::PARASAIL);
            if got == f[2] {
                same += 1;
            } else {
                if shown < 3 {
                    shown += 1;
                    eprintln!(
                        "differs on {}x{}:\n  sg    : {}\n  global: {}",
                        f[0].len(),
                        f[1].len(),
                        f[2],
                        got
                    );
                }
                differ += 1;
            }
            n += 1;
        }
        eprintln!("global vs parasail sg: {same} identical, {differ} different, of {n}");
    }
}
// Removed on the way across from isONcorrect: the whole `blockalign` test module.
//
// It measured how close block-aligner's SIMD banded alignment came to this exact
// path, and whether the difference survived `guard::fix_correction`. Both halves
// are isONcorrect-specific: block-aligner was on its default path for the
// structural-overcorrection guard, and isONform has neither the guard nor the
// dependency. The exact semi-global path below is what isONform needs, and its
// own oracle test against recorded `parasail.sg_trace_scan_16` output is what
// verifies it.

// Removed on the way across from isONcorrect: the `one_case` probe module, a
// one-off that read two sequences from files and compared semi-global against
// `global_affine`. Both the comparison target and the guard question behind it
// are gone.
