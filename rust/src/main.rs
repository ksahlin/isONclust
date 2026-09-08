//! isONclust -- Rust port. See PORTING.md.
//!
//! Only the CLI is ported so far; the clustering stages are not written yet, so
//! a valid invocation exits 3 saying so rather than silently producing nothing.

mod align;
mod blockalign;
mod cli;
mod cluster;
mod fastq;
mod minimizers;
mod p_emp;
mod packed;
mod parallelize;
mod parasail;
#[cfg(feature = "parasail-ffi")]
mod parasail_ffi;
mod phred;
mod pyfloat;
mod pyround;
mod simd;
mod sorting;
mod sweep;
mod text;
mod wfa;

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
        cli::Outcome::WriteFastq(wf) => write_fastq(&wf),
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
    // With no stage requested and --t 1, the port now runs end to end.
    if stage.is_empty() {
        return run_pipeline(&args, &outfolder);
    }

    // `minimizers` dumps the same format as `bench/dump_reference.py --stage
    // minimizers`, so the two can be diffed directly. It takes an ALREADY
    // SORTED fastq, because that is the order `reads_to_clusters` iterates and
    // the stage under test is the minimizer selection, not the sort.
    if stage == "minimizers" {
        return dump_minimizers(&args);
    }

    // `mapping` replays recorded get_best_cluster calls. The sweep is stateful
    // -- the minimizer database grows as reads become representatives -- so the
    // reference's calls are captured by wrapping the live driver
    // (bench/dump_reference.py --stage mapping) and replayed here.
    if stage == "mapping" {
        return replay_mapping(&args);
    }

    // `parasail` replays the alignments the reference actually performed.
    // isONform verified its parasail port against isONcorrect's parameters
    // (match 4, mismatch -8, open 12); isONclust uses match 2, mismatch -2 and
    // an opening penalty of 2..5, which can reach different tie-breaking paths.
    if stage == "parasail" {
        return replay_parasail();
    }

    // Compare candidate aligners against the exact one on real recorded calls.
    if stage == "aligners" {
        return compare_aligners(&args);
    }

    if stage != "sort" {
        eprintln!("isONclust: the clustering stages are not ported yet.");
        eprintln!("The CLI and the sorting stage are done; see PORTING.md 'Port status'.");
        eprintln!("Run the sorting stage alone with ISONCLUST_STAGE=sort.");
        return ExitCode::from(3);
    }

    ExitCode::from(match run_sort_stage(&args, &outfolder) {
        Ok(_) => 0,
        Err(code) => code,
    })
}

/// The sorting stage. Returns a process exit code; 0 on success.
/// Returns the number of reads written to `sorted.fastq`, so the clustering
/// stage can size its vector exactly instead of growing it by doubling.
/// `Ok(None)` means an existing sorted file was reused, so the count is unknown.
fn run_sort_stage(args: &cli::Args, outfolder: &str) -> Result<Option<usize>, u8> {
    let k = args.k as usize;
    let sorted_path = std::path::Path::new(&outfolder).join("sorted.fastq");
    let log_path = std::path::Path::new(&outfolder).join("logfile.txt");

    // The reference opens logfile.txt with mode 'w' as its FIRST action, before
    // deciding whether to do any work -- which is why --use_old_sorted_file
    // leaves it empty (PORTING.md, Finding 8).
    if let Err(e) = std::fs::write(&log_path, "") {
        eprintln!("isONclust: cannot write logfile: {e}");
        return Err(1);
    }
    if args.use_old_sorted_file && sorted_path.exists() {
        println!("Using already existing sorted file in specified directory, in not intended, specify different outfolder or delete the current file.");
        return Ok(None);
    }

    let path = match (&args.fastq, &args.flnc, &args.ccs) {
        (Some(f), _, _) => f.clone(),
        _ => {
            eprintln!("isONclust: the --ccs/--flnc BAM path is not ported (PORTING.md, Scope)");
            return Err(3);
        }
    };
    // Two streaming passes, keeping 32 bytes per read instead of its bases.
    //
    // Pass one scores every record and records where it sits in the input and how
    // long its output line will be. The vector is then sorted by score, which
    // fixes each surviving record's byte offset in `sorted.fastq`. Pass two
    // streams the input again and writes each record straight to its place with
    // `write_all_at`, so the input is read sequentially and only the output is
    // addressed out of order.
    //
    // Holding `acc`, `seq` and `qual` for every read was 1.73 GB on
    // SIRV_real_full -- and once the clustering stage's sequences were packed, it
    // was the largest thing in the whole run. See PORTING.md, "Memory: measured".
    struct SortRec {
        score: f64,
        error_rate: f64,
        /// Index of this record in the input, so pass two can find it again.
        ordinal: u32,
        /// Bytes `sorted_fastq_record` will write for it.
        out_len: u32,
    }
    let mut recs: Vec<SortRec> = Vec::new();
    // Finding 11: the reference crashes with `TypeError: 'NoneType' object is
    // not iterable` when a record has no quality -- which happens when the file
    // has no trailing newline. Reproduce the failure, with an explanation. The
    // reference reports the first such read, so keep the first and carry on.
    let mut no_qual: Option<String> = None;
    let mut ordinal: u32 = 0;
    if let Err(e) = fastq::for_each_file(std::path::Path::new(&path), |r| {
        let this = ordinal;
        ordinal += 1;
        if r.qual.is_none() && no_qual.is_none() {
            no_qual = Some(r.name.clone());
        }
        if let Some(sc) = sorting::score_record(&r, k, args.quality_threshold) {
            let qual = r.qual.as_deref().unwrap_or("");
            recs.push(SortRec {
                score: sc.score,
                error_rate: sc.error_rate,
                ordinal: this,
                out_len: sorting::sorted_fastq_record_len(&r.name, sc.score, &r.seq, qual) as u32,
            });
        }
    }) {
        // The reference dies with a traceback here; a message is better and
        // this path is not in the byte-identity contract.
        eprintln!("isONclust: cannot read {path}: {e}");
        return Err(1);
    }
    if let Some(name) = no_qual {
        eprintln!(
            "isONclust: read '{name}' has no quality values. The reference crashes here with"
        );
        eprintln!("TypeError: 'NoneType' object is not iterable. The usual cause is a fastq");
        eprintln!("with no trailing newline on its last line. See PORTING.md, Finding 11.");
        return Err(1);
    }

    let nr_scored = recs.len();
    sorting::sort_by_score(&mut recs, |r| r.score);

    // Where each surviving record lands in the output, indexed by input ordinal.
    // `u64::MAX` marks a record the quality filter dropped.
    let mut place: Vec<u64> = vec![u64::MAX; ordinal as usize];
    let mut total: u64 = 0;
    for r in &recs {
        place[r.ordinal as usize] = total;
        total += u64::from(r.out_len);
    }
    // The score has to be available again in pass two to rebuild the header, and
    // recomputing it would mean a second compensated sum over the quality string.
    let mut score_of_ordinal: Vec<f64> = vec![0.0; ordinal as usize];
    for r in &recs {
        score_of_ordinal[r.ordinal as usize] = r.score;
    }
    // `logfile_contents` sorts its input, so the order here does not matter.
    let mut rates: Vec<f64> = recs.iter().map(|r| r.error_rate).collect();
    drop(recs);

    {
        use std::os::unix::fs::FileExt;
        let f = match std::fs::File::create(&sorted_path) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("isONclust: cannot write sorted.fastq: {e}");
                return Err(1);
            }
        };
        if let Err(e) = f.set_len(total) {
            eprintln!("isONclust: cannot write sorted.fastq: {e}");
            return Err(1);
        }
        let mut ordinal: u32 = 0;
        let mut failed: Option<std::io::Error> = None;
        if let Err(e) = fastq::for_each_file(std::path::Path::new(&path), |r| {
            let this = ordinal as usize;
            ordinal += 1;
            let at = place[this];
            if at == u64::MAX || failed.is_some() {
                return;
            }
            let line = sorting::sorted_fastq_record(
                &r.name,
                score_of_ordinal[this],
                &r.seq,
                r.qual.as_deref().unwrap_or(""),
            );
            if let Err(e) = f.write_all_at(line.as_bytes(), at) {
                failed = Some(e);
            }
        }) {
            eprintln!("isONclust: cannot re-read {path}: {e}");
            return Err(1);
        }
        if let Some(e) = failed {
            eprintln!("isONclust: cannot write sorted.fastq: {e}");
            return Err(1);
        }
    }
    println!(
        "{} reads passed quality critera (avg phred Q val over {} and length > 2*k) and will be clustered.",
        nr_scored,
        pyfloat::repr(args.quality_threshold)
    );

    match sorting::logfile_contents(&mut rates) {
        Some(contents) => {
            if let Err(e) = std::fs::write(&log_path, contents) {
                eprintln!("isONclust: cannot write logfile: {e}");
                return Err(1);
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
            return Err(1);
        }
    }
    Ok(Some(nr_scored))
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

/// Replay recorded `get_best_cluster` calls and emit the same `RES` lines.
///
/// The dump path comes from `ISONCLUST_MAPPING_DUMP`. Format is
/// `bench/dump_reference.py --stage mapping`'s: a `CALL` line, one `CAND` line
/// per candidate in the reference's dict order, then `RES`. Only `RES` is
/// emitted here, so a diff against the reference's own `RES` lines is the check.
fn replay_mapping(args: &cli::Args) -> ExitCode {
    use std::collections::HashMap;

    let path = match std::env::var("ISONCLUST_MAPPING_DUMP") {
        Ok(p) => p,
        Err(_) => {
            eprintln!("isONclust: ISONCLUST_STAGE=mapping needs ISONCLUST_MAPPING_DUMP=<file>");
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
    let table = match p_emp::Table::select(args.k, args.w) {
        Some(t) => t,
        None => {
            // Finding 13: the reference builds an empty dict and dies with
            // KeyError on the first lookup.
            eprintln!(
                "isONclust: no empirical probabilities within +-2 of --w {} for --k {}.",
                args.w, args.k
            );
            return ExitCode::from(1);
        }
    };

    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::new(stdout.lock());

    let mut read_cl_id = 0usize;
    let mut seq_len = 0usize;
    let mut n_minimizers = 0usize;
    let mut hits = cluster::Hits::default();
    let mut reps: HashMap<usize, (String, f64)> = HashMap::new();

    let parse_list = |s: &str| -> Vec<usize> {
        if s.is_empty() {
            Vec::new()
        } else {
            s.split(',')
                .map(|x| x.parse().expect("integer list"))
                .collect()
        }
    };

    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        match f.first().copied() {
            Some("CALL") => {
                read_cl_id = f[1].parse().expect("read id");
                seq_len = f[2].parse().expect("seq len");
                n_minimizers = f[3].parse().expect("minimizer count");
                let e: f64 = f[4].parse().expect("error rate");
                hits = cluster::Hits::default();
                reps.clear();
                reps.insert(read_cl_id, (String::new(), e));
            }
            Some("CAND") => {
                let cl_id: usize = f[1].parse().expect("cluster id");
                let e: f64 = f[2].parse().expect("error rate");
                let indices = parse_list(f[3]);
                let positions = parse_list(f[4]);
                let acc = f.get(5).copied().unwrap_or("").to_string();
                hits.order.push(cl_id);
                hits.by_cluster
                    .insert(cl_id, cluster::HitList { indices, positions });
                reps.insert(cl_id, (acc, e));
            }
            Some("RES") => {
                struct ReplayReps<'a>(&'a HashMap<usize, (String, f64)>);
                impl cluster::Representatives for ReplayReps<'_> {
                    fn acc(&self, id: usize) -> &str {
                        &self.0[&id].0
                    }
                    fn error_rate(&self, id: usize) -> f64 {
                        self.0[&id].1
                    }
                }
                let r = cluster::get_best_cluster(
                    read_cl_id,
                    seq_len,
                    &hits,
                    n_minimizers,
                    &ReplayReps(&reps),
                    &table,
                    args.min_shared,
                    args.min_fraction,
                    args.min_prob_no_hits,
                    args.mapped_threshold,
                );
                let _ = writeln!(
                    out,
                    "RES\t{}\t{}\t{}",
                    r.best_cluster_id,
                    r.nr_shared_kmers,
                    pyfloat::repr(r.mapped_ratio)
                );
            }
            _ => {}
        }
    }
    let _ = out.flush();
    ExitCode::SUCCESS
}

/// Replay recorded `parasail_block_alignment` calls.
///
/// Dump path from `ISONCLUST_PARASAIL_DUMP`; format is
/// `bench/dump_reference.py --stage parasail`'s:
/// `PARA<TAB>open<TAB>k<TAB>match_id<TAB>s1<TAB>s2<TAB>cigar<TAB>ratio`.
/// Emits the same line with the port's cigar and ratio, so a plain diff is the
/// check -- both the alignment path chosen and the ratio computed from it.
fn replay_parasail() -> ExitCode {
    let path = match std::env::var("ISONCLUST_PARASAIL_DUMP") {
        Ok(p) => p,
        Err(_) => {
            eprintln!("isONclust: ISONCLUST_STAGE=parasail needs ISONCLUST_PARASAIL_DUMP=<file>");
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
    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::new(stdout.lock());
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.first().copied() != Some("PARA") || f.len() < 8 {
            continue;
        }
        let open: i32 = f[1].parse().expect("opening penalty");
        let k: usize = f[2].parse().expect("k");
        let match_id: i64 = f[3].parse().expect("match_id");
        let (s1, s2) = (f[4].as_bytes(), f[5].as_bytes());
        let aln = parasail::semiglobal(s1, s2, blockalign::scoring(open));
        let block = blockalign::parasail_block_alignment(s1, s2, k, match_id, open);
        let _ = writeln!(
            out,
            "PARA\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            open,
            k,
            match_id,
            f[4],
            f[5],
            aln.cigar,
            pyfloat::repr(block.alignment_ratio)
        );
    }
    let _ = out.flush();
    ExitCode::SUCCESS
}

/// Whether to run in exactly-the-reference mode.
///
/// **True by default, and that is the whole point of this port**: the
/// specification is byte-identity, so every approximation is opt-in. isONform's
/// equivalent defaults the other way, because its specification is accuracy --
/// the flag is carried across with `wfa.rs` and inverted here deliberately.
///
/// `ISONCLUST_FAITHFUL=0` turns the experimental aligners on. Nothing in the
/// tested contract runs with it set; it exists for `ISONCLUST_STAGE=aligners`.
pub fn faithful() -> bool {
    !matches!(
        std::env::var("ISONCLUST_FAITHFUL").as_deref(),
        Ok("0") | Ok("false")
    )
}

/// Strip the score suffix the sorting stage appended:
/// `"_".join(acc.split("_")[:-1])`.
pub fn strip_score(acc: &str) -> &str {
    match acc.rfind('_') {
        Some(i) => &acc[..i],
        // The reference's join of an empty list is "", not the original string.
        None => "",
    }
}

/// The score the sorting stage appended, recovered with `float(acc.split("_")[-1])`.
fn score_of(acc: &str) -> f64 {
    acc.rsplit('_')
        .next()
        .and_then(|s| s.parse().ok())
        .unwrap_or(f64::NAN)
}

/// The whole pipeline: sort, cluster (single core or batched), write.
fn run_pipeline(args: &cli::Args, outfolder: &str) -> ExitCode {
    let k = args.k as usize;
    let table = match p_emp::Table::select(args.k, args.w) {
        Some(t) => t,
        None => {
            // Finding 13: the reference builds an empty dict and then dies with
            // KeyError: (0.01, 0.01) partway through clustering.
            eprintln!(
                "isONclust: no empirical minimizer probabilities within +-2 of --w {} for --k {}.",
                args.w, args.k
            );
            eprintln!("The reference reaches this too, and fails with KeyError: (0.01, 0.01).");
            return ExitCode::from(1);
        }
    };

    // --- the sorting stage, as ISONCLUST_STAGE=sort does it ---
    let sorted_count = match run_sort_stage(args, outfolder) {
        Ok(n) => n,
        Err(code) => return ExitCode::from(code),
    };
    let sorted_path = std::path::Path::new(outfolder).join("sorted.fastq");
    // The reference re-reads sorted.fastq rather than reusing the in-memory
    // array, and the score is recovered from the accession. Reproduced, because
    // the accession the sweep sees is the one WITH the score suffix.
    //
    // Streamed straight into `SweepRead`. This used to slurp the file, parse it
    // into `Vec<Record>`, copy that into `Vec<Scored>` and copy *that* into
    // `Vec<SweepRead>` -- four resident copies of every base and quality score,
    // where one will do. The two intermediates were never read for anything
    // else. See PORTING.md, "Memory: measured".
    // Sized from the sort stage's own count, so the vector never reallocates:
    // growing 1.3M entries by doubling has the old and the new buffer live at
    // once at each step. `None` is the --use_old_sorted_file path, where no
    // count was produced.
    let mut sweep_reads: Vec<sweep::SweepRead> = match sorted_count {
        Some(n) => Vec::with_capacity(n),
        None => Vec::new(),
    };
    // Counted so a divergence cannot happen silently: 2-bit packing turns any
    // non-ACGT base into an `A`, which changes minimizer selection and the
    // alignment path. Every corpus here is pure ACGT, so this stays zero and the
    // harness cannot see the divergence -- Finding 5 again, hence the report.
    let mut substituted = 0usize;
    if let Err(e) = fastq::for_each_file(&sorted_path, |r| {
        let score = score_of(&r.name);
        let (seq, sub) = packed::PackedSeq::from_bytes(r.seq.as_bytes());
        substituted += sub;
        // Everything the clustering needs from the quality string, computed here
        // so the string itself never becomes resident. Both expressions are
        // written exactly as their consumers wrote them, so the f64s are
        // bit-identical to what the reference produces. See `SweepRead`.
        let qual = r.qual.unwrap_or_default();
        let qual_b = qual.as_bytes();
        let hp_error_rate = sweep::compressed_error_rate(r.seq.as_bytes(), qual_b);
        let err_per_base = blockalign::expected_errors(qual_b) / r.seq.len() as f64;
        sweep_reads.push(sweep::SweepRead {
            id: sweep_reads.len(),
            prev_batch_index: 0,
            acc: r.name.into(),
            seq,
            hp_error_rate,
            err_per_base,
            score,
        });
    }) {
        eprintln!("isONclust: cannot read {}: {e}", sorted_path.display());
        return ExitCode::from(1);
    }

    // Captured before the sweep consumes the reads; this is `reads.len()` as it
    // was when the two intermediate vectors still existed.
    if substituted > 0 {
        eprintln!(
            "isONclust: warning: {substituted} non-ACGT bases were read as 'A'. Sequences are\n\
             2-bit packed, which cannot represent a fifth symbol, so output for this input\n\
             will NOT match the Python reference. See rust/src/packed.rs."
        );
    }
    let nr_reads = sweep_reads.len();

    let params = sweep::SweepParams {
        k,
        w: args.w as usize,
        min_shared: args.min_shared,
        min_fraction: args.min_fraction,
        min_prob_no_hits: args.min_prob_no_hits,
        mapped_threshold: args.mapped_threshold,
        aligned_threshold: args.aligned_threshold,
    };
    let (clusters, representatives, stats) = if args.nr_cores > 1 {
        // --t > 1 is a DIFFERENT ALGORITHM, not a parallelised one. See
        // parallelize.rs and PORTING.md Finding 3.
        let bt = match parallelize::BatchType::parse(&args.batch_type) {
            Some(b) => b,
            None => {
                // Finding 15: the CLI documents "weighted", which no branch
                // implements; the reference produces no batches and dies with
                // "ValueError: Number of processes must be at least 1".
                eprintln!(
                    "isONclust: --batch_type {:?} is not implemented.",
                    args.batch_type
                );
                eprintln!("Use total_nt, nr_reads or read_lengths_squared. Note the help text's");
                eprintln!("\"weighted\" is not implemented in the reference either.");
                return ExitCode::from(1);
            }
        };
        let r = parallelize::parallel_clustering(
            &sweep_reads,
            args.nr_cores as usize,
            bt,
            &table,
            params,
            &sorted_path,
        );
        // The per-iteration files parallel mode writes.
        for (i, (pre, origins)) in r.intermediates.iter().enumerate() {
            let dir = std::path::Path::new(outfolder).join((i + 1).to_string());
            if let Err(e) = std::fs::create_dir_all(&dir) {
                eprintln!("isONclust: cannot create {}: {e}", dir.display());
                return ExitCode::from(1);
            }
            if let Err(e) = std::fs::write(dir.join("pre_clusters.csv"), pre)
                .and_then(|_| std::fs::write(dir.join("cluster_origins.csv"), origins))
            {
                eprintln!("isONclust: cannot write intermediates: {e}");
                return ExitCode::from(1);
            }
        }
        (
            r.clusters,
            r.representatives,
            (
                r.mapped_passed,
                r.aln_passed,
                r.aln_called,
                r.skipped_short,
                r.times,
            ),
        )
    } else {
        // Left empty: `reads_to_clusters` creates each read's own cluster and
        // representative when it first sees the read, and removes them again as
        // soon as the read maps. Pre-populating meant building a 1.3M-entry
        // table and 1.3M single-element Vecs, then discarding 99.6% of them.
        let clusters = sweep::OrderedClusters::default();
        let reps: rustc_hash::FxHashMap<usize, sweep::ReadInfo> = rustc_hash::FxHashMap::default();

        let res = sweep::reads_to_clusters(
            clusters,
            reps,
            &sweep_reads,
            cluster::MinimizerDatabase::new(),
            1,
            &table,
            params,
        );
        (
            res.clusters,
            res.representatives,
            (
                res.mapped_passed,
                res.aln_passed,
                res.aln_called,
                res.skipped_short,
                res.times,
            ),
        )
    };

    // --- write output, ordered by (cluster size, representative score) desc ---
    //
    // `sorted(..., reverse=True)` is stable and does NOT reverse ties, so equal
    // (size, score) pairs keep dict insertion order -- which after the
    // reassignment step is ascending cluster id. Hence the third key.
    let mut order: Vec<usize> = clusters.order.clone();
    order.sort_by(|a, b| {
        let ka = (clusters.map[a].len(), representatives[a].score);
        let kb = (clusters.map[b].len(), representatives[b].score);
        kb.0.cmp(&ka.0)
            .then(kb.1.partial_cmp(&ka.1).expect("scores are finite"))
            .then(a.cmp(b))
    });

    // The representatives' quality strings, which the clustering stage no longer
    // holds. One sequential pass over sorted.fastq for the survivors only.
    let want: rustc_hash::FxHashSet<usize> = order.iter().copied().collect();
    let rep_quals = match quals_for(&sorted_path, &want) {
        Ok(q) => q,
        Err(e) => {
            eprintln!("isONclust: cannot re-read {}: {e}", sorted_path.display());
            return ExitCode::from(1);
        }
    };

    let mut clusters_out = String::new();
    let mut origins_out = String::new();
    let mut nontrivial = 0usize;
    for (output_cl_id, c_id) in order.iter().enumerate() {
        let rep = &representatives[c_id];
        origins_out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            output_cl_id,
            strip_score(&rep.acc),
            String::from_utf8_lossy(&rep.seq.to_bytes()),
            rep_quals.get(c_id).map(String::as_str).unwrap_or(""),
            pyfloat::repr(rep.score),
            pyfloat::repr(rep.error_rate.unwrap_or(f64::NAN)),
        ));
        // Members are read ids; the accession and its score come from the
        // sorted read array. `SweepRead::score` was itself set by `score_of` on
        // the same accession at load, so this is the identical f64 the reference
        // parses back out, and the sort is unchanged. It is deliberately NOT
        // total -- ties keep member-list order, which is ascending read id.
        let mut members: Vec<u32> = clusters.map[c_id].clone();
        members.sort_by(|a, b| {
            sweep_reads[*b as usize]
                .score
                .partial_cmp(&sweep_reads[*a as usize].score)
                .expect("scores are finite")
        });
        for id in &members {
            clusters_out.push_str(&format!(
                "{}\t{}\n",
                output_cl_id,
                strip_score(&sweep_reads[*id as usize].acc)
            ));
        }
        if clusters.map[c_id].len() > 1 {
            nontrivial += 1;
        }
    }

    let cp = std::path::Path::new(outfolder).join("final_clusters.tsv");
    let op = std::path::Path::new(outfolder).join("final_cluster_origins.tsv");
    if let Err(e) = std::fs::write(&cp, clusters_out).and_then(|_| std::fs::write(&op, origins_out))
    {
        eprintln!("isONclust: cannot write output: {e}");
        return ExitCode::from(1);
    }
    let (mapped_passed, aln_passed, aln_called, skipped_short, times) = stats;
    if std::env::var("ISONCLUST_PROFILE").is_ok() {
        times.report(&format!("--t {}", args.nr_cores));
    }

    println!("Total number of reads iterated through:{}", nr_reads);
    println!("Passed mapping criteria:{}", mapped_passed);
    println!("Passed alignment criteria in this process:{}", aln_passed);
    println!(
        "Total calls to alignment mudule in this process:{}",
        aln_called
    );
    if skipped_short > 0 {
        // The reference prints one line per skipped read; a count says the same
        // thing without burying the summary. stdout is not in the contract.
        println!("skipped {skipped_short} reads whose homopolymer-compressed length was under k");
    }
    println!("Nr clusters larger than 1: {}", nontrivial);
    println!("Nr clusters (all): {}", order.len());
    ExitCode::SUCCESS
}

/// The `write_fastq` subcommand: split a clustering into per-cluster fastq files.
///
/// Two details that matter:
///
/// * **Cluster ids stay strings.** The reference reads them from the file and
///   uses them unparsed as the filename (`str(cl_id) + ".fastq"`), so a file
///   holding `007` would produce `007.fastq`. Parsing to an integer here would
///   quietly rename it.
/// * **The output order is the order ids first appear in the clusters file**,
///   because the reference iterates a `defaultdict`. It decides only which file
///   is created first, but it is free to preserve.
///
/// Accessions in the clusters file have had their score suffix stripped, and are
/// looked up in the ORIGINAL fastq -- not `sorted.fastq` -- so they match the
/// input's own read names.
/// Quality strings for a set of read ids, recovered by streaming `sorted.fastq`.
///
/// The clustering stage does not keep quality strings resident -- see
/// `sweep::SweepRead` -- and the two origins writers are the only things that
/// need them, for the surviving representatives only: 579 of 1 295 814 on
/// SIRV_real_full. A read's id is its ordinal in `sorted.fastq`, so one
/// sequential pass finds them all.
fn quals_for(
    sorted_path: &std::path::Path,
    want: &rustc_hash::FxHashSet<usize>,
) -> std::io::Result<rustc_hash::FxHashMap<usize, String>> {
    let mut out: rustc_hash::FxHashMap<usize, String> = rustc_hash::FxHashMap::default();
    let mut ordinal = 0usize;
    fastq::for_each_file(sorted_path, |r| {
        let this = ordinal;
        ordinal += 1;
        if want.contains(&this) {
            out.insert(this, r.qual.unwrap_or_default());
        }
    })?;
    Ok(out)
}

/// Read the record the index points at, and confirm it is the one asked for.
///
/// `Ok(None)` means the index had no entry, or the entry pointed at a different
/// read -- a hash collision. Either way the caller falls back to a scan. The
/// check is free: the record has to be parsed anyway to be re-emitted.
fn fetch_record(
    src: &std::fs::File,
    index: &rustc_hash::FxHashMap<u64, (u64, u32)>,
    key: u64,
    acc: &str,
    raw: &mut Vec<u8>,
) -> std::io::Result<Option<fastq::Record>> {
    use std::os::unix::fs::FileExt;
    let Some(&(at, n)) = index.get(&key) else {
        return Ok(None);
    };
    raw.resize(n as usize, 0);
    src.read_exact_at(raw, at)?;
    let text = String::from_utf8_lossy(raw);
    let mut rec = None;
    // Re-parsed with the same parser that produced the range, so there is no
    // second interpretation of the bytes to get out of step.
    fastq::for_each(text.split_inclusive('\n'), |r| rec = Some(r));
    Ok(rec.filter(|r| r.name == acc))
}

/// Find one record by accession by reading the whole file.
///
/// Only reached on a 64-bit hash collision, which for 100M reads has about a
/// 3e-4 chance of happening at all. Correct rather than fast.
fn scan_for(path: &std::path::Path, acc: &str) -> std::io::Result<Option<fastq::Record>> {
    let mut found = None;
    fastq::for_each_file(path, |r| {
        if found.is_none() && r.name == acc {
            found = Some(r);
        }
    })?;
    Ok(found)
}

/// One sequential scan of the fastq, writing every cluster as reads arrive.
///
/// Valid only when the fastq presents each cluster's reads in the order the
/// clusters file lists them, which holds for `sorted.fastq` because both are
/// score-descending with the same tie-break. `Ok(false)` means it did not hold
/// -- a read arrived out of turn, or a cluster finished short -- and the caller
/// should use the offset path. Nothing partial is left behind that the offset
/// path will not truncate.
///
/// The rank guard doubles as the collision check: two accessions sharing a hash
/// put one read in the wrong cluster, which shows up either as an out-of-turn
/// rank or as a short cluster at the end.
fn write_sequential(
    fastq_path: &std::path::Path,
    outfolder: &str,
    order: &[String],
    expected: &[u32],
    writable: &[bool],
    plan: &rustc_hash::FxHashMap<u64, (u32, u32)>,
) -> Result<bool, String> {
    use std::io::Write;

    let mut writers: Vec<Option<std::io::BufWriter<std::fs::File>>> =
        (0..order.len()).map(|_| None).collect();
    for (ci, name) in order.iter().enumerate() {
        if !writable[ci] {
            continue;
        }
        let path = std::path::Path::new(outfolder).join(format!("{name}.fastq"));
        match std::fs::File::create(&path) {
            Ok(f) => writers[ci] = Some(std::io::BufWriter::with_capacity(1 << 16, f)),
            // Out of descriptors despite the budget: let the caller seek instead.
            Err(e) if e.kind() == std::io::ErrorKind::Other || e.raw_os_error() == Some(24) => {
                return Ok(false)
            }
            Err(e) => return Err(format!("cannot write {}: {e}", path.display())),
        }
    }

    let mut next: Vec<u32> = vec![0; order.len()];
    let mut out_of_order = false;
    let mut io_err: Option<String> = None;
    fastq::for_each_file(fastq_path, |r| {
        if out_of_order || io_err.is_some() {
            return;
        }
        // `sorted.fastq` names carry the appended score that `final_clusters.tsv`
        // strips, so try the name as-is first -- which is what the original reads
        // file has -- and only then the stripped form. `strip_score` cannot be
        // applied unconditionally: it removes everything after the last `_`, and
        // real accessions end in things like `_strand=+`.
        let (key, name) = match plan.get(&hash_acc_for(&r.name)) {
            Some(v) => (*v, r.name.as_str()),
            None => {
                let stripped = strip_score(&r.name);
                match plan.get(&hash_acc_for(stripped)) {
                    Some(v) => (*v, stripped),
                    None => return,
                }
            }
        };
        let (ci, rank) = key;
        let ci = ci as usize;
        if !writable[ci] {
            return;
        }
        if rank != next[ci] {
            out_of_order = true;
            return;
        }
        let w = writers[ci].as_mut().expect("writable clusters are opened");
        if let Err(e) = write!(
            w,
            "@{}\n{}\n+\n{}\n",
            name,
            r.seq,
            r.qual.unwrap_or_default()
        ) {
            io_err = Some(format!("cannot write cluster fastq: {e}"));
            return;
        }
        next[ci] += 1;
    })
    .map_err(|e| format!("cannot read {}: {e}", fastq_path.display()))?;

    if let Some(e) = io_err {
        return Err(e);
    }
    if out_of_order {
        return Ok(false);
    }
    // Every cluster must have received exactly what the clusters file promised.
    for ci in 0..order.len() {
        if writable[ci] && next[ci] != expected[ci] {
            return Ok(false);
        }
    }
    for w in writers.iter_mut().flatten() {
        w.flush()
            .map_err(|e| format!("cannot write cluster fastq: {e}"))?;
    }
    Ok(true)
}

#[cfg(test)]
mod write_fastq_tests {
    /// `strip_score` removes everything after the last `_`, so it must not be
    /// applied to an accession that has no appended score -- real ONT names end
    /// in things like `_strand=+`. The sequential path relies on trying the raw
    /// name first for exactly this reason.
    #[test]
    fn strip_score_would_mangle_an_unscored_accession() {
        let scored = "read_45_abc/1_strand=+_779.1486123878182";
        assert_eq!(crate::strip_score(scored), "read_45_abc/1_strand=+");
        let unscored = "read_45_abc/1_strand=+";
        assert_eq!(crate::strip_score(unscored), "read_45_abc/1");
        assert_ne!(crate::strip_score(unscored), unscored);
    }
}

/// The accession hash both paths key on.
fn hash_acc_for(acc: &str) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h = std::collections::hash_map::DefaultHasher::new();
    acc.hash(&mut h);
    h.finish()
}

fn write_fastq(wf: &cli::WriteFastqArgs) -> ExitCode {
    let (clusters_path, fastq_path, outfolder) = match (&wf.clusters, &wf.fastq, &wf.outfolder) {
        (Some(c), Some(f), Some(o)) => (c, f, o),
        _ => {
            // The reference dies with a TypeError from os.path.join(None, ...).
            eprintln!(
                "isONclust write_fastq: --clusters, --fastq and --outfolder are all required"
            );
            return ExitCode::from(1);
        }
    };

    // Three streaming passes, and nothing per-read is held except a hash.
    //
    // The join here is between the clusters file, which names reads by
    // accession, and the fastq, which stores them in a different order, so
    // *something* has to map accession to record. Holding the accessions is what
    // this used to do and it does not scale: at 100M reads the clusters file's
    // strings plus a string-keyed index come to about 27 GB. Keying by a 64-bit
    // hash of the accession instead is 2.2 GB at that size.
    //
    // The hash makes this approximate, so it is checked: the record is parsed on
    // read-back anyway, and its name is compared with the accession asked for. A
    // mismatch, or an accession the index never saw, falls through to a scan of
    // the fastq, which is correct if slow. With 100M reads the chance of any
    // collision at all is about 3e-4.
    let hash_acc = hash_acc_for;

    // Pass 1: for every line, which cluster it belongs to and its rank inside
    // that cluster, keyed by a hash of the accession. Plus each cluster's name
    // and size. The accessions themselves are not kept: at 100M reads the
    // clusters file's strings alone are 11 GB.
    let mut idx_of: rustc_hash::FxHashMap<String, u32> = rustc_hash::FxHashMap::default();
    let mut order: Vec<String> = Vec::new();
    let mut expected: Vec<u32> = Vec::new();
    let mut plan: rustc_hash::FxHashMap<u64, (u32, u32)> = rustc_hash::FxHashMap::default();
    {
        use std::io::BufRead;
        let f = match std::fs::File::open(clusters_path) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("isONclust: cannot read {clusters_path}: {e}");
                return ExitCode::from(1);
            }
        };
        for line in std::io::BufReader::with_capacity(1 << 20, f).lines() {
            let line = match line {
                Ok(l) => l,
                Err(e) => {
                    eprintln!("isONclust: cannot read {clusters_path}: {e}");
                    return ExitCode::from(1);
                }
            };
            // `line.strip().split()` -- any whitespace, and blank lines vanish.
            let mut it = line.split_whitespace();
            let (cl_id, acc) = match (it.next(), it.next()) {
                (Some(a), Some(b)) => (a, b),
                _ => continue,
            };
            let ci = match idx_of.get(cl_id) {
                Some(&i) => i,
                None => {
                    let i = u32::try_from(order.len()).expect("fewer than 4G clusters");
                    order.push(cl_id.to_string());
                    expected.push(0);
                    idx_of.insert(cl_id.to_string(), i);
                    i
                }
            };
            let rank = expected[ci as usize];
            expected[ci as usize] += 1;
            plan.insert(hash_acc(acc), (ci, rank));
        }
    }
    let writable: Vec<bool> = expected.iter().map(|&c| (c as i64) >= wf.n).collect();
    let n_writable = writable.iter().filter(|w| **w).count();

    if let Err(e) = std::fs::create_dir_all(outfolder) {
        eprintln!("isONclust: cannot create {outfolder}: {e}");
        return ExitCode::from(1);
    }

    // The sequential path, when few enough clusters are being written that every
    // output file can be held open at once.
    //
    // `sorted.fastq` is score-descending, and each cluster's member list in
    // `final_clusters.tsv` is score-descending with the same tie-break, so a
    // cluster's reads appear in the fastq in exactly the order the clusters file
    // lists them -- checked on 178 multi-read clusters, no exceptions. When that
    // holds, one sequential scan can write every cluster with no seeking and no
    // byte-offset index, which matters when the fastq is far larger than page
    // cache: 377 GB at 100M PacBio reads, where the random-access path degrades
    // to real device I/O.
    //
    // It does not always hold: `--fastq` is documented as the *original* reads,
    // whose order is arbitrary. So the ranks recorded above are used as a guard.
    // A read arriving out of turn, or a cluster that ends up short, means the
    // input was not in cluster order -- or that two accessions collided on their
    // hash -- and the offset path runs instead. Either way the result is exact.
    const HANDLE_BUDGET: usize = 4096;
    if n_writable > 0 && n_writable <= HANDLE_BUDGET {
        match write_sequential(
            std::path::Path::new(fastq_path),
            outfolder,
            &order,
            &expected,
            &writable,
            &plan,
        ) {
            Ok(true) => {
                println!("Wrote clusters to separate fastq files.");
                return ExitCode::SUCCESS;
            }
            Ok(false) => { /* not in cluster order; fall through */ }
            Err(e) => {
                eprintln!("isONclust: {e}");
                return ExitCode::from(1);
            }
        }
    }
    // Not needed by the offset path, and holding both indexes at once would
    // double the footprint.
    drop(plan);

    // Pass 2: index the fastq. The byte ranges come from the parser, not a scan
    // for `@`, because a quality line can begin with `@`.
    let mut index: rustc_hash::FxHashMap<u64, (u64, u32)> = rustc_hash::FxHashMap::default();
    if let Err(e) = fastq::for_each_file_indexed(std::path::Path::new(fastq_path), |r, at, n| {
        index.insert(hash_acc(&r.name), (at, n));
    }) {
        eprintln!("isONclust: cannot read {fastq_path}: {e}");
        return ExitCode::from(1);
    }

    let src = match std::fs::File::open(fastq_path) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("isONclust: cannot read {fastq_path}: {e}");
            return ExitCode::from(1);
        }
    };
    // Pass 3: walk the clusters file again and write as we go. Records go
    // straight to the output rather than into a per-cluster String -- the
    // largest cluster on a 1.3M-read corpus holds 258 386 reads, which was
    // 172 MB of buffer.
    {
        use std::io::BufRead;
        use std::io::Write;
        let _ = &order;
        let f = match std::fs::File::open(clusters_path) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("isONclust: cannot read {clusters_path}: {e}");
                return ExitCode::from(1);
            }
        };
        // Which cluster files have been opened, so a clusters file that does not
        // group its lines appends instead of truncating.
        let mut opened: rustc_hash::FxHashSet<String> = rustc_hash::FxHashSet::default();
        let mut current: Option<(String, std::io::BufWriter<std::fs::File>)> = None;
        let mut raw: Vec<u8> = Vec::new();

        for line in std::io::BufReader::with_capacity(1 << 20, f).lines() {
            let line = match line {
                Ok(l) => l,
                Err(e) => {
                    eprintln!("isONclust: cannot read {clusters_path}: {e}");
                    return ExitCode::from(1);
                }
            };
            let mut it = line.split_whitespace();
            let (cl_id, acc) = match (it.next(), it.next()) {
                (Some(a), Some(b)) => (a, b),
                _ => continue,
            };
            let Some(&ci) = idx_of.get(cl_id) else {
                continue;
            };
            if !writable[ci as usize] {
                continue;
            }
            if current.as_ref().map(|(id, _)| id.as_str()) != Some(cl_id) {
                if let Some((_, mut w)) = current.take() {
                    if let Err(e) = w.flush() {
                        eprintln!("isONclust: cannot write cluster fastq: {e}");
                        return ExitCode::from(1);
                    }
                }
                let path = std::path::Path::new(outfolder).join(format!("{cl_id}.fastq"));
                let append = !opened.insert(cl_id.to_string());
                let file = std::fs::OpenOptions::new()
                    .write(true)
                    .create(true)
                    .append(append)
                    .truncate(!append)
                    .open(&path);
                let file = match file {
                    Ok(f) => f,
                    Err(e) => {
                        eprintln!("isONclust: cannot write {}: {e}", path.display());
                        return ExitCode::from(1);
                    }
                };
                current = Some((
                    cl_id.to_string(),
                    std::io::BufWriter::with_capacity(1 << 20, file),
                ));
            }

            let rec = match fetch_record(&src, &index, hash_acc(acc), acc, &mut raw) {
                Ok(Some(r)) => r,
                Ok(None) => match scan_for(std::path::Path::new(fastq_path), acc) {
                    Ok(Some(r)) => r,
                    Ok(None) => {
                        // The reference raises KeyError here.
                        eprintln!(
                            "isONclust: read {acc:?} is in {clusters_path} but not in {fastq_path}"
                        );
                        return ExitCode::from(1);
                    }
                    Err(e) => {
                        eprintln!("isONclust: cannot read {fastq_path}: {e}");
                        return ExitCode::from(1);
                    }
                },
                Err(e) => {
                    eprintln!("isONclust: cannot read {fastq_path}: {e}");
                    return ExitCode::from(1);
                }
            };
            let (_, w) = current.as_mut().expect("a writer is open");
            if let Err(e) = write!(
                w,
                "@{}\n{}\n+\n{}\n",
                acc,
                rec.seq,
                rec.qual.unwrap_or_default()
            ) {
                eprintln!("isONclust: cannot write cluster fastq: {e}");
                return ExitCode::from(1);
            }
        }
        if let Some((_, mut w)) = current.take() {
            if let Err(e) = w.flush() {
                eprintln!("isONclust: cannot write cluster fastq: {e}");
                return ExitCode::from(1);
            }
        }
    }
    println!("Wrote clusters to separate fastq files.");
    ExitCode::SUCCESS
}

/// Replay the recorded alignments through every candidate aligner and report
/// what each costs and what each changes.
///
/// The metric that matters is **not** CIGAR equality. isONclust reads one number
/// out of the alignment -- the fraction of windows with enough matches -- and
/// compares it against `--aligned_threshold`. Two different optimal paths can
/// give the same verdict, and a tiny ratio difference on the wrong side of the
/// threshold changes a read's cluster. So this reports, in order of increasing
/// relevance: time, CIGAR agreement, ratio agreement, and **verdict agreement**.
///
/// That ordering is isONform's finding 41 -- gate an aligner swap on verdicts,
/// not on scores -- applied to this tool's actual decision.
fn compare_aligners(args: &cli::Args) -> ExitCode {
    let path = match std::env::var("ISONCLUST_PARASAIL_DUMP") {
        Ok(p) => p,
        Err(_) => {
            eprintln!("isONclust: ISONCLUST_STAGE=aligners needs ISONCLUST_PARASAIL_DUMP=<file>");
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

    struct Case {
        open: i32,
        k: usize,
        match_id: i64,
        s1: Vec<u8>,
        s2: Vec<u8>,
        ref_cigar: String,
        ref_ratio: f64,
    }
    let mut cases: Vec<Case> = Vec::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.first().copied() != Some("PARA") || f.len() < 8 {
            continue;
        }
        cases.push(Case {
            open: f[1].parse().expect("open"),
            k: f[2].parse().expect("k"),
            match_id: f[3].parse().expect("match_id"),
            s1: f[4].as_bytes().to_vec(),
            s2: f[5].as_bytes().to_vec(),
            ref_cigar: f[6].to_string(),
            ref_ratio: f[7].parse().expect("ratio"),
        });
    }
    if cases.is_empty() {
        eprintln!("isONclust: no PARA records in {path}");
        return ExitCode::from(1);
    }

    let threshold = args.aligned_threshold;
    println!("# {} recorded alignments from {}", cases.len(), path);
    println!("# verdict = (alignment_ratio >= --aligned_threshold {threshold})");
    println!(
        "{:<14} {:>9} {:>7} {:>12} {:>12} {:>12} {:>9}",
        "ALIGNER", "SECS", "REL", "CIGAR_SAME", "RATIO_SAME", "VERDICT_SAME", "DECLINED"
    );

    // The ratio computation, shared, so only the aligner differs.
    let ratio_from =
        |ops: &[align::CigarOp], s1: &[u8], s2: &[u8], k: usize, match_id: i64| -> f64 {
            let (a1, a2) = match align::ops_to_seq(ops, s1, s2) {
                Some(x) => x,
                None => return f64::NAN,
            };
            let matches: Vec<u8> = a1
                .iter()
                .zip(a2.iter())
                .map(|(x, y)| u8::from(x == y))
                .collect();
            let head = matches.len().min(k);
            let mut current: i64 = matches[..head].iter().map(|x| i64::from(*x)).sum();
            let mut aligned: i64 = i64::from(current >= match_id);
            for (leaving, &new_state) in matches.iter().zip(matches.iter().skip(k)) {
                current = current - i64::from(*leaving) + i64::from(new_state);
                aligned += i64::from(current >= match_id);
            }
            aligned as f64 / s1.len() as f64
        };

    let mut baseline_secs = 0f64;
    for which in ["parasail(exact)", "block-aligner", "wfa2"] {
        let t0 = std::time::Instant::now();
        let (mut cig_same, mut ratio_same, mut verdict_same, mut declined) = (0usize, 0, 0, 0);
        for c in &cases {
            let sc = blockalign::scoring(c.open);
            let (cigar, ops): (String, Vec<align::CigarOp>) = match which {
                "parasail(exact)" => {
                    let a = parasail::semiglobal(&c.s1, &c.s2, sc);
                    (a.cigar, a.ops)
                }
                "block-aligner" => {
                    let mut ops = Vec::new();
                    match simd::semiglobal_ops(&c.s1, &c.s2, sc, &mut ops) {
                        Some(_) => (align::encode_cigar(&ops), ops),
                        None => {
                            declined += 1;
                            let a = parasail::semiglobal(&c.s1, &c.s2, sc);
                            (a.cigar, a.ops)
                        }
                    }
                }
                _ => match wfa::semiglobal(&c.s1, &c.s2, sc) {
                    Some(a) => (a.cigar, a.ops),
                    None => {
                        declined += 1;
                        let a = parasail::semiglobal(&c.s1, &c.s2, sc);
                        (a.cigar, a.ops)
                    }
                },
            };
            let ratio = ratio_from(&ops, &c.s1, &c.s2, c.k, c.match_id);
            if cigar == c.ref_cigar {
                cig_same += 1;
            }
            if ratio == c.ref_ratio {
                ratio_same += 1;
            }
            if (ratio >= threshold) == (c.ref_ratio >= threshold) {
                verdict_same += 1;
            }
        }
        let secs = t0.elapsed().as_secs_f64();
        if which == "parasail(exact)" {
            baseline_secs = secs;
        }
        let n = cases.len() as f64;
        println!(
            "{:<14} {:>9.2} {:>6.2}x {:>11.4}% {:>11.4}% {:>11.4}% {:>9}",
            which,
            secs,
            baseline_secs / secs,
            100.0 * cig_same as f64 / n,
            100.0 * ratio_same as f64 / n,
            100.0 * verdict_same as f64 / n,
            declined
        );
    }
    ExitCode::SUCCESS
}
