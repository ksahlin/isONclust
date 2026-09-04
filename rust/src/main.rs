//! isONclust -- Rust port. See PORTING.md.
//!
//! Only the CLI is ported so far; the clustering stages are not written yet, so
//! a valid invocation exits 3 saying so rather than silently producing nothing.

mod cli;
mod text;

use std::io::Write;
use std::process::ExitCode;

/// Flags dropped from the port. PORTING.md, *Scope*: the consensus feature is
/// out of scope entirely. Refusing by name matters -- a pipeline that passes
/// `--consensus`, gets exit 0 and no `consensus_references.fasta` is worse off
/// than one that fails loudly, and a bare "unrecognised argument" is not
/// actionable.
const DROPPED: &[&str] = &[
    "--consensus",
    "--medaka",
    "--abundance_ratio",
    "--rc_identity_threshold",
];

fn main() -> ExitCode {
    let argv: Vec<String> = std::env::args().skip(1).collect();

    // Checked before parsing, so the message is about the dropped flag rather
    // than about whatever else is wrong with the command line.
    for a in &argv {
        let name = a.split('=').next().unwrap_or(a);
        if let Some(d) = DROPPED.iter().find(|d| **d == name) {
            eprintln!("isONclust: error: {} is not supported in the Rust port.", d);
            eprintln!(
                "The consensus feature (--consensus, --abundance_ratio, --rc_identity_threshold, --medaka)"
            );
            eprintln!(
                "was dropped: isONclust clusters reads, and consensus of a cluster is isONcorrect's job."
            );
            return ExitCode::from(2);
        }
    }

    match cli::parse(&argv) {
        cli::Outcome::Exit0(out) => {
            print!("{}", out);
            let _ = std::io::stdout().flush();
            ExitCode::SUCCESS
        }
        cli::Outcome::UsageError(msg) => {
            eprint!("{}", text::USAGE);
            eprintln!("isONclust: error: {}", msg);
            ExitCode::from(2)
        }
        cli::Outcome::Message(msg, code) => {
            print!("{}", msg);
            let _ = std::io::stdout().flush();
            ExitCode::from(code as u8)
        }
        cli::Outcome::Run(_args) => {
            eprintln!("isONclust: the clustering stages are not ported yet.");
            eprintln!("The CLI is complete and verified; see PORTING.md 'Port status'.");
            ExitCode::from(3)
        }
        cli::Outcome::WriteFastq(_wf) => {
            eprintln!("isONclust: write_fastq is not ported yet.");
            ExitCode::from(3)
        }
    }
}
