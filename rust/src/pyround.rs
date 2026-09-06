//! Python's `round(x, 2)`.
//!
//! `p_shared_minimizer_empirical` rounds both error rates to two decimals and
//! uses the result as a dict key, so this decides which of the 15 probability
//! buckets a read lands in -- and therefore which reads cluster together.
//!
//! Python's `round(x, n)` is **correctly rounded on the exact binary value**,
//! with ties (which need an exactly representable midpoint) going to even. It is
//! not "multiply by 100, round, divide": that double-rounds and disagrees.
//! Measured on the reference:
//!
//! | x | `round(x, 2)` | why |
//! |---|---|---|
//! | 0.005 | **0.01** | the double is 0.005000000000000000104…, just above the midpoint |
//! | 0.015 | **0.01** | the double is 0.01499999999999999944…, just below |
//! | 0.025 | **0.03** | just above |
//! | 0.045 | **0.04** | just below |
//! | 0.125 | **0.12** | exactly representable, so the tie breaks to even |
//!
//! Rust's `{:.2}` formatting is also exact and also breaks ties to even, so
//! formatting and parsing back reproduces it. `rust/tests/pyround_oracle.rs`
//! checks that against CPython over a large sample rather than trusting it.

/// `round(x, 2)`.
pub fn round2(x: f64) -> f64 {
    if !x.is_finite() {
        return x;
    }
    // `{:.2}` renders the exact decimal expansion of the double, rounded to two
    // places, ties to even -- the same rule Python applies.
    format!("{:.2}", x)
        .parse()
        .expect("a formatted float always parses")
}

#[cfg(test)]
mod tests {
    use super::round2;

    /// Values taken from the reference, not from reasoning about them.
    #[test]
    fn matches_cpython_on_the_awkward_cases() {
        assert_eq!(round2(0.005), 0.01);
        assert_eq!(round2(0.015), 0.01);
        assert_eq!(round2(0.025), 0.03);
        assert_eq!(round2(0.035), 0.04);
        assert_eq!(round2(0.045), 0.04);
        assert_eq!(round2(0.055), 0.06);
        assert_eq!(round2(0.065), 0.07);
        assert_eq!(round2(0.075), 0.07);
        assert_eq!(round2(0.085), 0.09);
        assert_eq!(round2(0.095), 0.1);
        assert_eq!(round2(0.125), 0.12);
        assert_eq!(round2(0.135), 0.14);
    }

    /// The naive form disagrees, which is why this module exists.
    #[test]
    fn the_obvious_implementation_is_wrong() {
        let naive = |x: f64| (x * 100.0).round() / 100.0;
        assert_ne!(naive(0.015), round2(0.015));
        assert_ne!(naive(0.045), round2(0.045));
    }

    #[test]
    fn ordinary_values_round_normally() {
        assert_eq!(round2(0.1234), 0.12);
        assert_eq!(round2(0.0), 0.0);
        assert_eq!(round2(0.19), 0.19);
        assert_eq!(round2(1.0), 1.0);
    }
}
