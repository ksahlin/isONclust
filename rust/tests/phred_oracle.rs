//! The frozen phred tables must still match the live reference.
//!
//! `src/phred.rs` is generated from the reference environment because Rust's
//! `powf` and CPython's `**` disagree by one ULP on `%` (phred 4). A frozen
//! table is only correct while the environment that produced it is the one in
//! use, so this checks -- rather than trusting the comment.

use std::io::Write;
use std::process::Command;

// The module is compiled standalone here, so items this test does not
// exercise look dead. They are not -- they are used by the binary.
#[allow(dead_code)]
#[path = "../src/phred.rs"]
mod phred;

#[test]
fn frozen_tables_match_the_reference() {
    if std::env::var("ISONCLUST_SKIP_PYTHON_ORACLE").is_ok() {
        eprintln!("SKIPPED by ISONCLUST_SKIP_PYTHON_ORACLE -- phred tables unverified here");
        return;
    }
    let py = std::env::var("REF_PYTHON").unwrap_or_else(|_| {
        format!(
            "{}/miniforge3/envs/isonclust-ref/bin/python",
            std::env::var("HOME").unwrap_or_default()
        )
    });

    let mut input = String::new();
    for (i, cap, unc) in phred::all() {
        input.push_str(&format!(
            "{}\t{:016x}\t{:016x}\n",
            i,
            cap.to_bits(),
            unc.to_bits()
        ));
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
    assert!(
        out.status.success() && stdout.contains("OK"),
        "frozen phred tables no longer match the reference:\n{stdout}\n\
         Regenerate src/phred.rs -- and note this means the reference environment changed, \
         so every golden needs re-recording too."
    );
    eprintln!("{}", stdout.trim());
}

const PY_CHECK: &str = r#"
import struct, sys
bad = 0
n = 0
for line in sys.stdin:
    i_s, cap_hex, unc_hex = line.rstrip("\n").split("\t")
    i = int(i_s)
    want_unc = 10 ** (-(i - 33) / 10.0)
    want_cap = min(want_unc, 0.79433)
    got_cap = struct.unpack(">d", bytes.fromhex(cap_hex))[0]
    got_unc = struct.unpack(">d", bytes.fromhex(unc_hex))[0]
    n += 2
    for label, want, got in (("capped", want_cap, got_cap), ("uncapped", want_unc, got_unc)):
        if struct.pack(">d", want) != struct.pack(">d", got):
            bad += 1
            if bad <= 8:
                print(f"MISMATCH char={i} ({chr(i)!r}) {label}: python={want!r} frozen={got!r}")
if bad:
    print(f"{bad} of {n} entries differ")
    sys.exit(1)
print(f"OK {n} frozen phred entries identical to the reference")
"#;
