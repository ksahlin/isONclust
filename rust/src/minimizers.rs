//! `cluster.get_kmer_minimizers`, ported.
//!
//! Minimizers are chosen by **lexicographic order on the k-mer string**, not by
//! a hash. Three properties of the reference decide bytes here, and two of them
//! look like bugs:
//!
//! 1. **Ties go to the FIRST position in the window** (`list(window).index`).
//!    isONcorrect uses `rindex` and takes the last; do not carry that habit
//!    across.
//! 2. **The window can run past the end of the sequence.** It is built as
//!    `[seq[i:i+k] for i in range(w - k + 1)]` with no length check, and Python
//!    slicing silently yields short strings, then empty ones. The only guard
//!    upstream is `len(hpol) < k`, so any read whose compressed length falls in
//!    `[k, w)` reaches this. An empty string is lexicographically smallest, so
//!    it wins its window -- and an empty minimizer is a database key that
//!    matches every other short read. Live on real data: 445 of 19 972
//!    `droso_20k` reads at the default `--k 15 --w 50`. See PORTING.md,
//!    Finding 4.
//! 3. **The same `(minimizer, position)` can be emitted twice.** When the k-mer
//!    leaving the window merely *equals* the current minimum by value, the
//!    minimum is recomputed and appended again, which can land on the same
//!    position. `get_all_hits` counts per minimizer index, so duplicates inflate
//!    hit counts.
//!
//! All three are reproduced. Verified against the reference by
//! `bench/equivalence.sh stage minimizers`.

/// A k-mer slice, clamped to the end of the sequence exactly as Python's
/// `seq[i:i+k]` clamps -- and empty once `i` is past the end.
#[inline]
fn kmer(seq: &[u8], i: usize, k: usize) -> &[u8] {
    if i >= seq.len() {
        return &[];
    }
    let end = (i + k).min(seq.len());
    &seq[i..end]
}

/// `get_kmer_minimizers(seq, k_size, w_size)`.
///
/// Returns `(minimizer, position)` pairs in the reference's order. Positions are
/// offsets into `seq`, which the caller has already homopolymer-compressed.
pub fn get_kmer_minimizers<'a>(seq: &'a [u8], k: usize, w_size: usize) -> Vec<(&'a [u8], usize)> {
    // `w = w_size - k_size`. The CLI guarantees w_size >= k, so this cannot
    // wrap, but saturating_sub keeps the function total for direct callers.
    let w = w_size.saturating_sub(k);

    // The reference builds this without checking the length, so entries beyond
    // the end are short or empty. Reproduced.
    let mut window: std::collections::VecDeque<&[u8]> = (0..=w).map(|i| kmer(seq, i, k)).collect();

    let mut curr_min: &[u8] = *window
        .iter()
        .min()
        .expect("the window always holds at least one entry");
    let first_idx = window
        .iter()
        .position(|x| *x == curr_min)
        .expect("the minimum is in the window");
    let mut minimizers: Vec<(&[u8], usize)> = vec![(curr_min, first_idx)];

    // `range(w + 1, len(seq) - k_size + 1)`, which is empty when the sequence is
    // shorter than the first window.
    let upper = (seq.len() + 1).saturating_sub(k);
    for i in (w + 1)..upper {
        let new_kmer = kmer(seq, i, k);
        let discarded = window.pop_front().expect("window is non-empty");
        window.push_back(new_kmer);

        if discarded == curr_min {
            // The previous minimum has left the window *by value*, so recompute
            // brute force. Note this fires on equality, not identity, so a
            // repeated k-mer can re-emit the same pair.
            curr_min = *window.iter().min().expect("window is non-empty");
            let idx = window
                .iter()
                .position(|x| *x == curr_min)
                .expect("the minimum is in the window");
            minimizers.push((curr_min, idx + i - w));
        } else if new_kmer < curr_min {
            curr_min = new_kmer;
            minimizers.push((curr_min, i));
        }
    }
    minimizers
}

/// Homopolymer-compress, then take minimizers -- what `reads_to_clusters` does.
/// Returns `None` for a read the reference skips.
pub fn minimizers_for_read(
    seq: &[u8],
    k: usize,
    w: usize,
) -> Option<(Vec<u8>, Vec<(usize, usize)>)> {
    let hpol = crate::sorting::homopolymer_compress(seq);
    if hpol.len() < k {
        return None;
    }
    let ms = get_kmer_minimizers(&hpol, k, w);
    // Return offsets rather than borrowed slices so the compressed sequence can
    // be handed back alongside them.
    let spans = ms
        .iter()
        .map(|(m, pos)| (*pos, m.len()))
        .collect::<Vec<_>>();
    Some((hpol, spans))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mins(seq: &str, k: usize, w: usize) -> Vec<(String, usize)> {
        get_kmer_minimizers(seq.as_bytes(), k, w)
            .into_iter()
            .map(|(m, p)| (String::from_utf8_lossy(m).into_owned(), p))
            .collect()
    }

    /// Finding 4, measured against the reference on an 18 nt sequence.
    #[test]
    fn the_window_runs_past_the_end_and_emits_empty_and_sub_k_minimizers() {
        let s = "ACGTACGTACGTACGTAC"; // 18 nt, no homopolymers
        assert_eq!(mins(s, 15, 50), vec![(String::new(), 18)]);
        assert_eq!(mins(s, 15, 20), vec![("ACGTACGTACGTAC".to_string(), 4)]);
        assert_eq!(mins(s, 13, 20), vec![("ACGTACGTACGTA".to_string(), 0)]);
    }

    /// An empty minimizer is a database key that matches every other short read.
    #[test]
    fn an_empty_minimizer_is_the_only_one_emitted() {
        let m = mins("ACGTACGTACGTACGTAC", 15, 50);
        assert_eq!(m.len(), 1);
        assert!(m[0].0.is_empty());
    }

    /// Ties resolve to the FIRST occurrence, not the last.
    #[test]
    fn ties_go_to_the_first_position() {
        // "AAAA" repeated: every 2-mer in the first window is "AA"
        let m = mins("AAAAAAAA", 2, 4);
        assert_eq!(m[0], ("AA".to_string(), 0));
    }

    #[test]
    fn a_window_of_one_makes_every_kmer_a_minimizer() {
        // w == k means w - k == 0, so the window holds a single k-mer
        let m = mins("ACGTACGT", 4, 4);
        let positions: Vec<usize> = m.iter().map(|x| x.1).collect();
        assert_eq!(positions, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn minimizers_are_non_decreasing_in_position() {
        let s = "ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA";
        for (k, w) in [(4usize, 8usize), (5, 10), (7, 15)] {
            let m = mins(s, k, w);
            let mut prev = 0usize;
            for (_, p) in &m {
                assert!(*p >= prev, "positions went backwards at {p} (k={k} w={w})");
                prev = *p;
            }
        }
    }

    #[test]
    fn short_reads_are_skipped_by_the_caller() {
        assert!(minimizers_for_read(b"AAAAAAAA", 15, 50).is_none());
        assert!(minimizers_for_read(b"ACGTACGTACGTACGTAC", 15, 50).is_some());
    }
}
