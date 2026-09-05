//! `pyround::round2` against CPython's `round(x, 2)`.
//!
//! The unit tests assert twelve values the author looked up. This asserts
//! agreement over a hundred thousand the author did not, concentrated on the
//! range the tool actually produces -- error rates in `[0, 0.3]` -- and on exact
//! hundredth-midpoints, where the rule is subtlest.

use std::io::Write;
use std::process::Command;

// Compiled standalone here, so items the test does not use look dead.
#[allow(dead_code)]
#[path = "../src/pyround.rs"]
mod pyround;

#[test]
fn round2_matches_cpython() {
    if std::env::var("ISONCLUST_SKIP_PYTHON_ORACLE").is_ok() {
        eprintln!("SKIPPED by ISONCLUST_SKIP_PYTHON_ORACLE -- pyround is unverified here");
        return;
    }
    let py = std::env::var("REF_PYTHON").unwrap_or_else(|_| {
        format!(
            "{}/miniforge3/envs/isonclust-ref/bin/python",
            std::env::var("HOME").unwrap_or_default()
        )
    });

    let mut vals: Vec<f64> = Vec::new();
    // every exact midpoint in range, where ties bite
    for h in 0..=300 {
        vals.push(h as f64 / 100.0 + 0.005);
        vals.push(h as f64 / 100.0);
    }
    // and a dense sweep of realistic error rates
    let mut x = 0x9E3779B97F4A7C15u64;
    for _ in 0..100_000 {
        x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
        vals.push((x >> 11) as f64 / 2f64.powi(53) * 0.3);
    }

    let mut input = String::with_capacity(vals.len() * 40);
    for v in &vals {
        input.push_str(&format!(
            "{:016x}\t{:016x}\n",
            v.to_bits(),
            pyround::round2(*v).to_bits()
        ));
    }

    let mut child = Command::new(&py)
        .arg("-c")
        .arg(PY_CHECK)
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
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
    assert!(out.status.success() && stdout.contains("OK"), "{stdout}");
    eprintln!("{}", stdout.trim());
}

const PY_CHECK: &str = r#"
import struct, sys
bad = n = 0
for line in sys.stdin:
    xh, rh = line.rstrip("\n").split("\t")
    x = struct.unpack(">d", bytes.fromhex(xh))[0]
    ours = struct.unpack(">d", bytes.fromhex(rh))[0]
    theirs = round(x, 2)
    n += 1
    if struct.pack(">d", ours) != struct.pack(">d", theirs):
        bad += 1
        if bad <= 10:
            print(f"MISMATCH x={x!r} python={theirs!r} rust={ours!r}")
if bad:
    print(f"{bad} of {n} differ"); sys.exit(1)
print(f"OK {n} values identical to CPython round(x, 2)")
"#;
