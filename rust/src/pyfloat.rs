//! Python's `str(float)` / `repr(float)`, reproduced exactly.
//!
//! This is not cosmetic. `get_sorted_fastq_for_cluster.main` writes
//! `acc + "_{0}".format(score)` into every read header of `sorted.fastq`, and
//! every later stage reads the score back with `float(acc.split("_")[-1])`. So
//! the formatted score is simultaneously an output byte and an input value, and
//! Rust's own `{}` does not agree with Python's:
//!
//! | value  | Python   | Rust `{}`          |
//! |--------|----------|--------------------|
//! | 1234.0 | `1234.0` | `1234`             |
//! | 1e-5   | `1e-05`  | `0.00001`          |
//! | 1e16   | `1e+16`  | `10000000000000000`|
//!
//! Both produce the same *shortest round-tripping digits*; they disagree only on
//! how those digits are laid out. So this takes Rust's digits and re-lays them
//! out Python's way.
//!
//! Python (CPython `format_float_short`, mode `'r'`) uses fixed notation when
//! `-4 < decpt <= 16` and exponential otherwise, where `decpt` is the position
//! of the decimal point relative to the digit string: `value = 0.d1d2..dn *
//! 10^decpt`. Measured boundaries: `1e15` -> `1000000000000000.0`, `1e16` ->
//! `1e+16`, `1e-4` -> `0.0001`, `1e-5` -> `1e-05`.

/// Format `v` exactly as Python's `str()` would.
pub fn repr(v: f64) -> String {
    if v.is_nan() {
        return "nan".to_string();
    }
    if v.is_infinite() {
        return if v > 0.0 { "inf" } else { "-inf" }.to_string();
    }

    let neg = v.is_sign_negative();
    let a = v.abs();

    if a == 0.0 {
        // Python keeps the sign of negative zero: str(-0.0) == '-0.0'
        return if neg { "-0.0" } else { "0.0" }.to_string();
    }

    // Rust's LowerExp gives the shortest round-tripping digits in the form
    // d[.ddd]e[-]dd, which is exactly the (digits, exponent) pair we need.
    let sci = format!("{:e}", a);
    let (mantissa, exp_str) = sci.split_once('e').expect("LowerExp always emits 'e'");
    let exp: i32 = exp_str
        .parse()
        .expect("LowerExp always emits an integer exponent");
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();

    // value = 0.d1d2..dn * 10^decpt
    let decpt = exp + 1;

    let body = if decpt <= -4 || decpt > 16 {
        format_exponential(&digits, exp)
    } else {
        format_fixed(&digits, decpt)
    };

    if neg {
        format!("-{}", body)
    } else {
        body
    }
}

fn format_exponential(digits: &str, exp: i32) -> String {
    let mut out = String::with_capacity(digits.len() + 6);
    out.push_str(&digits[..1]);
    if digits.len() > 1 {
        out.push('.');
        out.push_str(&digits[1..]);
    }
    out.push('e');
    // Python always writes the sign and pads the exponent to at least two
    // digits: 1e+16, 1e-05, 5e-324, 1.7976931348623157e+308.
    if exp < 0 {
        out.push('-');
    } else {
        out.push('+');
    }
    let mag = exp.unsigned_abs();
    if mag < 10 {
        out.push('0');
    }
    out.push_str(&mag.to_string());
    out
}

fn format_fixed(digits: &str, decpt: i32) -> String {
    let n = digits.len() as i32;
    if decpt <= 0 {
        // 0.0001 -> digits "1", decpt -3 -> "0." + "000" + "1"
        let mut out = String::from("0.");
        for _ in 0..(-decpt) {
            out.push('0');
        }
        out.push_str(digits);
        out
    } else if decpt >= n {
        // 1e15 -> digits "1", decpt 16 -> "1" + fifteen zeros + ".0"
        let mut out = String::from(digits);
        for _ in 0..(decpt - n) {
            out.push('0');
        }
        out.push_str(".0");
        out
    } else {
        let k = decpt as usize;
        format!("{}.{}", &digits[..k], &digits[k..])
    }
}

#[cfg(test)]
mod tests {
    use super::repr;

    /// The cases that motivated the module, taken from the reference's own
    /// behaviour rather than from memory.
    #[test]
    fn matches_python_on_the_documented_cases() {
        assert_eq!(repr(1234.0), "1234.0");
        assert_eq!(repr(1234.5678901234567), "1234.5678901234567");
        assert_eq!(repr(1e15), "1000000000000000.0");
        assert_eq!(repr(1e16), "1e+16");
        assert_eq!(repr(1e17), "1e+17");
        assert_eq!(repr(1e-4), "0.0001");
        assert_eq!(repr(1e-5), "1e-05");
        assert_eq!(repr(0.1 + 0.2), "0.30000000000000004");
        assert_eq!(repr(2409.0), "2409.0");
        assert_eq!(repr(1.0), "1.0");
        assert_eq!(repr(123456789012345.6), "123456789012345.6");
        assert_eq!(repr(1e22), "1e+22");
    }

    #[test]
    fn boundaries_are_where_python_puts_them() {
        // upper: fixed through 1e15, exponential from 1e16
        assert_eq!(repr(1e14), "100000000000000.0");
        assert_eq!(repr(1e15), "1000000000000000.0");
        assert_eq!(repr(1e16), "1e+16");
        // lower: fixed through 1e-4, exponential from 1e-5
        assert_eq!(repr(1e-3), "0.001");
        assert_eq!(repr(1e-4), "0.0001");
        assert_eq!(repr(1e-5), "1e-05");
        assert_eq!(repr(1e-6), "1e-06");
    }

    #[test]
    fn signs_zeros_and_specials() {
        assert_eq!(repr(0.0), "0.0");
        assert_eq!(repr(-0.0), "-0.0");
        assert_eq!(repr(-1234.0), "-1234.0");
        assert_eq!(repr(-1e-5), "-1e-05");
        assert_eq!(repr(f64::INFINITY), "inf");
        assert_eq!(repr(f64::NEG_INFINITY), "-inf");
        assert_eq!(repr(f64::NAN), "nan");
    }

    #[test]
    fn extremes_round_trip() {
        assert_eq!(repr(5e-324), "5e-324");
        assert_eq!(repr(f64::MAX), "1.7976931348623157e+308");
        assert_eq!(repr(f64::MIN_POSITIVE), "2.2250738585072014e-308");
    }

    /// Whatever this prints must parse back to the same bits, for every value.
    /// That is the property the reference actually depends on, because the
    /// score is written into the accession and parsed out again downstream.
    #[test]
    fn round_trips_through_parse() {
        let mut x = 1u64;
        for _ in 0..20000 {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let v = f64::from_bits(x);
            if !v.is_finite() {
                continue;
            }
            let s = repr(v);
            let back: f64 = s.parse().expect("must parse back");
            assert_eq!(back.to_bits(), v.to_bits(), "round trip failed for {}", s);
        }
    }
}
