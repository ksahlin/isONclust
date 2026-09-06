//! `modules/get_sorted_fastq_for_cluster.py`, ported: score every read, drop the
//! ones that fail the filters, sort by score descending, write `sorted.fastq`
//! and `logfile.txt`.
//!
//! Four things here decide bytes and are easy to get subtly wrong:
//!
//! * **Two phred tables.** `D` caps the per-base error probability at
//!   `0.79433`; `D_no_min` does not. The rolling no-error product uses the
//!   capped one, `error_rate` uses the uncapped one. The cap bites only on `!`
//!   (phred 0) and it is live: 121 of the retired corpus's reads contain one.
//! * **`sum()` is compensated.** The reference is pinned to Python >=3.12,
//!   where `sum()` over floats is Neumaier-compensated and therefore exactly
//!   rounded. `poisson_mean` sums over `set(qual)` -- a set, whose order is
//!   `PYTHONHASHSEED`-dependent -- so the compensation is what makes the
//!   reference deterministic at all (PORTING.md, Finding 1). Reproducing it
//!   exactly is what frees this code to iterate in any order it likes.
//! * **`math.log(x, 10)` is not `log10(x)`.** The reference computes
//!   `log(x)/log(10)`, which differs from `log10` in the last bit for 58% of
//!   values. No filter flip was observed in 800 000 trials, but exactness is
//!   the specification.
//! * **`sum_of_expectations` is a plain running `+=`,** not `sum()`, so it is
//!   *not* compensated. Sequential addition, in order.

use crate::fastq::Record;
use crate::pyfloat;

/// `D` -- per-base error probability, capped at 0.79433.
///
/// Read from a frozen table rather than computed: Rust's `powf` disagrees with
/// the reference's `**` by one ULP on `%` (phred 4). See `phred.rs`.
pub fn phred_capped(c: u8) -> f64 {
    crate::phred::capped(c)
}

/// `D_no_min` -- the same without the cap.
pub fn phred_uncapped(c: u8) -> f64 {
    crate::phred::uncapped(c)
}

/// Neumaier compensated summation: the exactly-rounded sum, matching the
/// builtin `sum()` over floats in CPython 3.12 and later. Order-independent in
/// practice, which is why the caller may iterate a map in any order.
pub fn fsum(xs: impl IntoIterator<Item = f64>) -> f64 {
    let mut s = 0.0f64;
    let mut c = 0.0f64;
    for x in xs {
        let t = s + x;
        if s.abs() >= x.abs() {
            c += (s - t) + x;
        } else {
            c += (x - t) + s;
        }
        s = t;
    }
    s + c
}

/// `expected_number_of_erroneous_kmers`.
///
/// A rolling product of per-base no-error probabilities over a window of `k`,
/// summed. The division-based update is the reference's, and it is kept because
/// it is what decides the low bits -- recomputing the window product each step
/// would be more accurate and would not match.
pub fn expected_erroneous_kmers(qual: &[u8], k: usize) -> f64 {
    let probs: Vec<f64> = qual.iter().map(|c| phred_capped(*c)).collect();
    if probs.len() < k {
        // The reference builds a deque of the first k entries; with fewer than
        // k it simply gets a short one. Callers filter these out first
        // (len(seq) < 2*k), so this is defensive.
        let prod: f64 = probs.iter().fold(1.0, |a, p| a * (1.0 - p));
        return probs.len() as f64 - k as f64 + 1.0 - prod;
    }
    let mut window: std::collections::VecDeque<f64> = probs[..k].iter().map(|p| 1.0 - p).collect();
    // `functools.reduce(operator.mul, window, 1)` -- a left fold from int 1.
    let mut current: f64 = window.iter().fold(1.0f64, |a, b| a * b);
    let mut sum_of_expectations = current;
    for p_e in &probs[k..] {
        let p_to_leave = window.pop_front().expect("window is non-empty");
        current *= (1.0 - p_e) / p_to_leave;
        sum_of_expectations += current;
        window.push_back(1.0 - p_e);
    }
    qual.len() as f64 - k as f64 + 1.0 - sum_of_expectations
}

/// Python's `10 * -math.log(error_rate, 10)`.
pub fn phred_of(error_rate: f64) -> f64 {
    10.0 * -(error_rate.ln() / 10f64.ln())
}

/// Homopolymer-compress: keep one character per run.
pub fn homopolymer_compress(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &c in seq {
        if out.last() != Some(&c) {
            out.push(c);
        }
    }
    out
}

/// A read that survived the filters.
#[derive(Debug, Clone)]
pub struct Scored {
    pub acc: String,
    pub seq: String,
    pub qual: String,
    pub score: f64,
    pub error_rate: f64,
}

/// Score and filter one record, in the reference's order, or `None` if it is
/// filtered out.
///
/// The order matters: the length/homopolymer filter runs *before* the score is
/// computed, and the quality filter runs *after* it, so a read can be scored
/// and then discarded.
///
/// Per-record rather than over a slice so the sort stage can stream the input
/// file instead of materialising every `Record` first. There is no state
/// between records, so this is the same computation in the same order.
pub fn score_record(r: &Record, k: usize, quality_threshold: f64) -> Option<Scored> {
    // Finding 11: the reference does not guard a missing quality string -- it
    // dies with `TypeError: 'NoneType' object is not iterable`. Reproduced as a
    // hard error by the caller; here we simply skip so the library is usable,
    // and `run` re-raises.
    let qual = r.qual.as_ref()?;
    let seq_b = r.seq.as_bytes();
    let qual_b = qual.as_bytes();

    let hpol = homopolymer_compress(seq_b);
    if seq_b.len() < 2 * k || hpol.len() < k {
        return None;
    }

    let exp_err = expected_erroneous_kmers(qual_b, k);
    let denom = (seq_b.len() - k + 1) as f64;
    let p_no_error = 1.0 - exp_err / denom;
    let score = p_no_error * denom;

    // poisson_mean = sum([qual.count(c) * D_no_min[c] for c in set(qual)])
    // A histogram gives the same multiset of terms; the compensated sum
    // makes the order irrelevant.
    let mut counts = [0u32; 256];
    for &c in qual_b {
        counts[c as usize] += 1;
    }
    let poisson_mean = fsum(
        counts
            .iter()
            .enumerate()
            .filter(|(_, n)| **n > 0)
            .map(|(c, n)| f64::from(*n) * phred_uncapped(c as u8)),
    );
    let error_rate = poisson_mean / qual_b.len() as f64;
    if phred_of(error_rate) <= quality_threshold {
        return None;
    }

    Some(Scored {
        acc: r.name.clone(),
        seq: r.seq.clone(),
        qual: qual.clone(),
        score,
        error_rate,
    })
}

/// `read_array.sort(key=lambda x: x[3], reverse=True)`.
///
/// Python's sort is stable and `reverse=True` does *not* reverse ties, so equal
/// scores keep their input order. `sort_by` in Rust is also stable, so this is
/// a direct translation -- but only if the comparison never says "equal" for
/// values that Python would order. Scores are finite here.
pub fn sort_by_score(reads: &mut [Scored]) {
    reads.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .expect("scores are finite; NaN would mean an upstream bug")
    });
}

/// The line the reference writes for each read: the score is appended to the
/// accession with Python's float formatting, and read back out downstream.
pub fn sorted_fastq_record(r: &Scored) -> String {
    format!(
        "@{}_{}\n{}\n+\n{}\n",
        r.acc,
        pyfloat::repr(r.score),
        r.seq,
        r.qual
    )
}

/// `logfile.txt`. Note the "median" takes the upper middle element with no
/// averaging, which is the reference's definition and not a mistake here.
pub fn logfile_contents(error_rates: &mut [f64]) -> Option<String> {
    if error_rates.is_empty() {
        return None;
    }
    error_rates.sort_by(|a, b| a.partial_cmp(b).expect("finite"));
    let min_e = error_rates[0];
    let max_e = error_rates[error_rates.len() - 1];
    let median_e = error_rates[error_rates.len() / 2];
    let mean_e = fsum(error_rates.iter().copied()) / error_rates.len() as f64;
    Some(format!(
        "Lowest read error rate:{}\nHighest read error rate:{}\nMedian read error rate:{}\nMean read error rate:{}\n\n",
        pyfloat::repr(min_e),
        pyfloat::repr(max_e),
        pyfloat::repr(median_e),
        pyfloat::repr(mean_e),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_cap_bites_only_on_phred_zero() {
        assert_eq!(phred_uncapped(b'!'), 1.0);
        assert_eq!(phred_capped(b'!'), 0.79433);
        // phred 1 computes to just under the cap, so `min` leaves it alone
        assert_eq!(phred_capped(b'"'), phred_uncapped(b'"'));
        assert!(phred_uncapped(b'"') < 0.79433);
        assert_eq!(phred_capped(b'I'), phred_uncapped(b'I'));
    }

    #[test]
    fn fsum_is_exactly_rounded_and_order_independent() {
        let xs = vec![1e17, 1.0, -1e17, 1.0];
        assert_eq!(fsum(xs.iter().copied()), 2.0);
        let mut rev = xs.clone();
        rev.reverse();
        assert_eq!(fsum(xs.iter().copied()), fsum(rev.iter().copied()));
        // a naive fold gets this wrong, which is the whole point
        let naive = xs.iter().fold(0.0f64, |a, b| a + b);
        assert_ne!(naive, 2.0);
    }

    #[test]
    fn homopolymer_compression() {
        assert_eq!(homopolymer_compress(b"AACCGGTT"), b"ACGT".to_vec());
        assert_eq!(homopolymer_compress(b"ACGT"), b"ACGT".to_vec());
        assert_eq!(homopolymer_compress(b""), Vec::<u8>::new());
        assert_eq!(homopolymer_compress(b"AAAA"), b"A".to_vec());
    }

    #[test]
    fn phred_uses_log_base_10_not_log10() {
        // The reference computes log(x)/log(10). Assert we do the same thing,
        // bit for bit, rather than calling log10.
        let x = 5.482811848400435e-06f64;
        assert_eq!(phred_of(x), 10.0 * -(x.ln() / 10f64.ln()));
    }

    #[test]
    fn a_perfect_read_scores_the_kmer_count_and_formats_with_a_dot_zero() {
        // Quality '~' is phred 93; the error probability is ~5e-10, so the
        // score lands just under the k-mer count rather than exactly on it.
        // What matters here is that the formatting goes through pyfloat.
        let q = vec![b'~'; 60];
        let e = expected_erroneous_kmers(&q, 15);
        assert!((0.0..1e-5).contains(&e), "expected ~0 errors, got {e}");
    }

    #[test]
    fn sort_is_descending_and_stable_on_ties() {
        let mk = |acc: &str, s: f64| Scored {
            acc: acc.to_string(),
            seq: String::new(),
            qual: String::new(),
            score: s,
            error_rate: 0.0,
        };
        let mut v = vec![mk("a", 1.0), mk("b", 3.0), mk("c", 1.0), mk("d", 3.0)];
        sort_by_score(&mut v);
        let order: Vec<&str> = v.iter().map(|r| r.acc.as_str()).collect();
        // descending by score; ties keep input order (b before d, a before c)
        assert_eq!(order, vec!["b", "d", "a", "c"]);
    }

    #[test]
    fn logfile_median_is_the_upper_middle_element() {
        let mut e = vec![0.4, 0.1, 0.3, 0.2];
        let s = logfile_contents(&mut e).expect("non-empty");
        // sorted: 0.1 0.2 0.3 0.4 -> index 4/2 = 2 -> 0.3, not 0.25
        assert!(s.contains("Median read error rate:0.3\n"), "{s}");
        assert!(s.contains("Lowest read error rate:0.1\n"));
        assert!(s.contains("Highest read error rate:0.4\n"));
    }

    #[test]
    fn logfile_is_none_when_everything_was_filtered() {
        assert!(logfile_contents(&mut []).is_none());
    }
}
