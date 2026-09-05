//! isONclust -- Rust port. See PORTING.md.
//!
//! Only the CLI is ported so far; the clustering stages are not written yet, so
//! a valid invocation exits 3 saying so rather than silently producing nothing.

mod cli;
mod fastq;
mod minimizers;
mod phred;
mod pyfloat;
mod sorting;
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
        cli::Outcome::Run(args) => run(*args),
        cli::Outcome::WriteFastq(_wf) => {
            eprintln!("isONclust: write_fastq is not ported yet.");
            ExitCode::from(3)
        }
    }
}

/// How far the port can currently go.
///
/// `ISONCLUST_STAGE=sort` runs the sorting stage and stops, so it can be
/// verified against the reference before the clustering exists. It is an
/// internal switch, deliberately an environment variable rather than a flag:
/// the CLI is a byte-for-byte contract and must not grow options the reference
/// does not have.
fn run(args: cli::Args) -> ExitCode {
    let outfolder = match &args.outfolder {
        Some(o) => o.clone(),
        None => {
            // The reference only creates the folder when --outfolder is given,
            // then fails on os.path.join(None, ...). Not reached by any golden.
            eprintln!("isONclust: --outfolder is required");
            return ExitCode::from(1);
        }
    };
    if let Err(e) = std::fs::create_dir_all(&outfolder) {
        eprintln!("isONclust: cannot create {outfolder}: {e}");
        return ExitCode::from(1);
    }

    let stage = std::env::var("ISONCLUST_STAGE").unwrap_or_default();

    // `minimizers` dumps the same format as `bench/dump_reference.py --stage
    // minimizers`, so the two can be diffed directly. It takes an ALREADY
    // SORTED fastq, because that is the order `reads_to_clusters` iterates and
    // the stage under test is the minimizer selection, not the sort.
    if stage == "minimizers" {
        return dump_minimizers(&args);
    }

    if stage != "sort" {
        eprintln!("isONclust: the clustering stages are not ported yet.");
        eprintln!("The CLI and the sorting stage are done; see PORTING.md 'Port status'.");
        eprintln!("Run the sorting stage alone with ISONCLUST_STAGE=sort.");
        return ExitCode::from(3);
    }

    let k = args.k as usize;
    let sorted_path = std::path::Path::new(&outfolder).join("sorted.fastq");
    let log_path = std::path::Path::new(&outfolder).join("logfile.txt");

    // The reference opens logfile.txt with mode 'w' as its FIRST action, before
    // deciding whether to do any work -- which is why --use_old_sorted_file
    // leaves it empty (PORTING.md, Finding 8).
    if let Err(e) = std::fs::write(&log_path, "") {
        eprintln!("isONclust: cannot write logfile: {e}");
        return ExitCode::from(1);
    }
    if args.use_old_sorted_file && sorted_path.exists() {
        println!("Using already existing sorted file in specified directory, in not intended, specify different outfolder or delete the current file.");
        return ExitCode::SUCCESS;
    }

    let path = match (&args.fastq, &args.flnc, &args.ccs) {
        (Some(f), _, _) => f.clone(),
        _ => {
            eprintln!("isONclust: the --ccs/--flnc BAM path is not ported (PORTING.md, Scope)");
            return ExitCode::from(3);
        }
    };
    let text = match std::fs::read_to_string(&path) {
        Ok(t) => t,
        Err(e) => {
            // The reference dies with a traceback here; a message is better and
            // this path is not in the byte-identity contract.
            eprintln!("isONclust: cannot read {path}: {e}");
            return ExitCode::from(1);
        }
    };

    let records = fastq::read(&text);
    // Finding 11: the reference crashes with `TypeError: 'NoneType' object is
    // not iterable` when a record has no quality -- which happens when the file
    // has no trailing newline. Reproduce the failure, with an explanation.
    if let Some(bad) = records.iter().find(|r| r.qual.is_none()) {
        eprintln!(
            "isONclust: read '{}' has no quality values. The reference crashes here with",
            bad.name
        );
        eprintln!("TypeError: 'NoneType' object is not iterable. The usual cause is a fastq");
        eprintln!("with no trailing newline on its last line. See PORTING.md, Finding 11.");
        return ExitCode::from(1);
    }

    let mut scored = sorting::score_reads(&records, k, args.quality_threshold);
    sorting::sort_by_score(&mut scored);

    let mut body = String::new();
    for r in &scored {
        body.push_str(&sorting::sorted_fastq_record(r));
    }
    if let Err(e) = std::fs::write(&sorted_path, &body) {
        eprintln!("isONclust: cannot write sorted.fastq: {e}");
        return ExitCode::from(1);
    }
    println!(
        "{} reads passed quality critera (avg phred Q val over {} and length > 2*k) and will be clustered.",
        scored.len(),
        pyfloat::repr(args.quality_threshold)
    );

    let mut rates: Vec<f64> = scored.iter().map(|r| r.error_rate).collect();
    match sorting::logfile_contents(&mut rates) {
        Some(contents) => {
            if let Err(e) = std::fs::write(&log_path, contents) {
                eprintln!("isONclust: cannot write logfile: {e}");
                return ExitCode::from(1);
            }
        }
        None => {
            // Matches the fix committed to the Python on master.
            let msg = format!(
                "No reads passed the quality filter (--q {}).\n",
                pyfloat::repr(args.quality_threshold)
            );
            let _ = std::fs::write(&log_path, &msg);
            eprintln!(
                "Error: no reads passed the quality filter (--q {}).",
                pyfloat::repr(args.quality_threshold)
            );
            eprintln!("Lower --q, or check that the input has quality values.");
            return ExitCode::from(1);
        }
    }
    ExitCode::SUCCESS
}

/// Dump `(read index, position, minimizer)` for every read, in file order.
///
/// Matches `bench/dump_reference.py --stage minimizers` byte for byte. A read
/// the reference skips (homopolymer-compressed length < k) emits one line with
/// position -1 and an empty minimizer, so the skip is part of the comparison
/// rather than an absence.
fn dump_minimizers(args: &cli::Args) -> ExitCode {
    let path = match &args.fastq {
        Some(f) => f.clone(),
        None => {
            eprintln!("isONclust: ISONCLUST_STAGE=minimizers needs --fastq <sorted.fastq>");
            return ExitCode::from(1);
        }
    };
    let text = match std::fs::read_to_string(&path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("isONclust: cannot read {path}: {e}");
            return ExitCode::from(1);
        }
    };
    let k = args.k as usize;
    let w = args.w as usize;

    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::new(stdout.lock());
    for (idx, rec) in fastq::read(&text).iter().enumerate() {
        match minimizers::minimizers_for_read(rec.seq.as_bytes(), k, w) {
            None => {
                let _ = writeln!(out, "{idx}\t-1\t");
            }
            Some((hpol, spans)) => {
                for (pos, len) in spans {
                    // The minimizer can be shorter than k, or empty, where the
                    // window ran past the end -- Finding 4.
                    let end = (pos + len).min(hpol.len());
                    let m = if pos >= hpol.len() {
                        &hpol[0..0]
                    } else {
                        &hpol[pos..end]
                    };
                    let _ = writeln!(out, "{idx}\t{pos}\t{}", String::from_utf8_lossy(m));
                }
            }
        }
    }
    let _ = out.flush();
    ExitCode::SUCCESS
}
