//! Differential oracle for `pyfloat::repr` against the live Python reference.
//!
//! The unit tests in `src/pyfloat.rs` assert what the author *believed* Python
//! does. This asserts what Python actually does, over values the author would
//! not have thought of. That distinction is the whole method: every hypothesis
//! stated before it was measured in the two prior ports turned out wrong.
//!
//! Requires an interpreter. If you genuinely have none, set
//! `ISONCLUST_SKIP_PYTHON_ORACLE=1` -- the skip is then explicit and visible in
//! the output, rather than the test quietly passing on a machine that never
//! ran it.

use std::io::Write;
use std::process::Command;

/// The same LCG as the unit test, so the two cover the same ground.
fn values(n: usize) -> Vec<f64> {
    let mut out = Vec::with_capacity(n);
    let mut x = 0x2545F4914F6CDD1Du64;
    // Values the tool actually produces: a score is
    // (1 - errors/(len-k+1)) * (len-k+1), so O(10^2..10^5) with a fraction.
    for i in 0..n {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let v = match i % 4 {
            // realistic scores
            0 => (x >> 11) as f64 / 2f64.powi(53) * 100000.0,
            // integral scores -- the case Rust's `{}` gets wrong, and which a
            // perfect read reaches (zero expected errors => score == len-k+1)
            1 => ((x >> 40) as f64).trunc(),
            // arbitrary bit patterns, including subnormals and huge values
            2 => f64::from_bits(x),
            // small magnitudes, to cross the 1e-5 boundary
            _ => (x >> 11) as f64 / 2f64.powi(53) * 1e-4,
        };
        if v.is_finite() {
            out.push(v);
        }
    }
    out
}

#[test]
fn repr_matches_cpython() {
    if std::env::var("ISONCLUST_SKIP_PYTHON_ORACLE").is_ok() {
        eprintln!("SKIPPED by ISONCLUST_SKIP_PYTHON_ORACLE -- pyfloat is unverified here");
        return;
    }
    let py = std::env::var("REF_PYTHON").unwrap_or_else(|_| {
        format!(
            "{}/miniforge3/envs/isonclust-ref/bin/python",
            std::env::var("HOME").unwrap_or_default()
        )
    });

    let vals = values(40000);

    // Send bits + our rendering; let Python check both directions.
    let mut input = String::with_capacity(vals.len() * 40);
    for v in &vals {
        input.push_str(&format!("{:016x}\t{}\n", v.to_bits(), isonclust_repr(*v)));
    }

    let mut child = Command::new(&py)
        .arg("-c")
        .arg(PY_CHECK)
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .unwrap_or_else(|e| {
            panic!(
                "could not run the reference interpreter at {py}: {e}\n\
                 Build it with bench/setup_reference_env.sh, point REF_PYTHON at it, \
                 or set ISONCLUST_SKIP_PYTHON_ORACLE=1 to skip explicitly."
            )
        });
    child
        .stdin
        .as_mut()
        .expect("stdin")
        .write_all(input.as_bytes())
        .expect("write");
    let out = child.wait_with_output().expect("wait");
    let stdout = String::from_utf8_lossy(&out.stdout);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        out.status.success(),
        "python oracle failed:\nstdout:\n{stdout}\nstderr:\n{stderr}"
    );
    assert!(
        stdout.contains("OK"),
        "mismatches against CPython:\n{stdout}"
    );
    eprintln!("{}", stdout.trim());
}

/// Re-exported through the binary crate's module path.
fn isonclust_repr(v: f64) -> String {
    // `src/pyfloat.rs` belongs to the binary crate, so it is included here
    // directly rather than imported. Keeping one copy of the source means the
    // oracle cannot drift from the implementation.
    include_pyfloat::repr(v)
}

// The module is compiled standalone here, so items this test does not
// exercise look dead. They are not -- they are used by the binary.
#[allow(dead_code)]
#[path = "../src/pyfloat.rs"]
mod include_pyfloat;

const PY_CHECK: &str = r#"
import struct, sys
bad = 0
n = 0
for line in sys.stdin:
    bits_hex, ours = line.rstrip("\n").split("\t")
    v = struct.unpack(">d", bytes.fromhex(bits_hex))[0]
    theirs = str(v)
    n += 1
    if ours != theirs:
        bad += 1
        if bad <= 10:
            print(f"MISMATCH bits={bits_hex} python={theirs!r} rust={ours!r}")
if bad:
    print(f"{bad} of {n} differ")
    sys.exit(1)
print(f"OK {n} values identical to CPython str(float)")
"#;
