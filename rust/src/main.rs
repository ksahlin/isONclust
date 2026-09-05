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
mod parallelize;
mod parasail;
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

    ExitCode::from(run_sort_stage(&args, &outfolder))
}

/// The sorting stage. Returns a process exit code; 0 on success.
fn run_sort_stage(args: &cli::Args, outfolder: &str) -> u8 {
    let k = args.k as usize;
    let sorted_path = std::path::Path::new(&outfolder).join("sorted.fastq");
    let log_path = std::path::Path::new(&outfolder).join("logfile.txt");

    // The reference opens logfile.txt with mode 'w' as its FIRST action, before
    // deciding whether to do any work -- which is why --use_old_sorted_file
    // leaves it empty (PORTING.md, Finding 8).
    if let Err(e) = std::fs::write(&log_path, "") {
        eprintln!("isONclust: cannot write logfile: {e}");
        return 1;
    }
    if args.use_old_sorted_file && sorted_path.exists() {
        println!("Using already existing sorted file in specified directory, in not intended, specify different outfolder or delete the current file.");
        return 0;
    }

    let path = match (&args.fastq, &args.flnc, &args.ccs) {
        (Some(f), _, _) => f.clone(),
        _ => {
            eprintln!("isONclust: the --ccs/--flnc BAM path is not ported (PORTING.md, Scope)");
            return 3;
        }
    };
    let text = match std::fs::read_to_string(&path) {
        Ok(t) => t,
        Err(e) => {
            // The reference dies with a traceback here; a message is better and
            // this path is not in the byte-identity contract.
            eprintln!("isONclust: cannot read {path}: {e}");
            return 1;
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
        return 1;
    }

    let mut scored = sorting::score_reads(&records, k, args.quality_threshold);
    sorting::sort_by_score(&mut scored);

    let mut body = String::new();
    for r in &scored {
        body.push_str(&sorting::sorted_fastq_record(r));
    }
    if let Err(e) = std::fs::write(&sorted_path, &body) {
        eprintln!("isONclust: cannot write sorted.fastq: {e}");
        return 1;
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
                return 1;
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
            return 1;
        }
    }
    0
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
    let rc = run_sort_stage(args, outfolder);
    if rc != 0 {
        return ExitCode::from(rc);
    }
    let sorted_path = std::path::Path::new(outfolder).join("sorted.fastq");
    let text = match std::fs::read_to_string(&sorted_path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("isONclust: cannot read {}: {e}", sorted_path.display());
            return ExitCode::from(1);
        }
    };
    // The reference re-reads sorted.fastq rather than reusing the in-memory
    // array, and the score is recovered from the accession. Reproduced, because
    // the accession the sweep sees is the one WITH the score suffix.
    let records = fastq::read(&text);
    let reads: Vec<sorting::Scored> = records
        .iter()
        .map(|r| sorting::Scored {
            acc: r.name.clone(),
            seq: r.seq.clone(),
            qual: r.qual.clone().unwrap_or_default(),
            score: score_of(&r.name),
            error_rate: f64::NAN,
        })
        .collect();

    let params = sweep::SweepParams {
        k,
        w: args.w as usize,
        min_shared: args.min_shared,
        min_fraction: args.min_fraction,
        min_prob_no_hits: args.min_prob_no_hits,
        mapped_threshold: args.mapped_threshold,
        aligned_threshold: args.aligned_threshold,
    };
    let sweep_reads: Vec<sweep::SweepRead> = reads
        .iter()
        .enumerate()
        .map(|(i, r)| sweep::SweepRead {
            id: i,
            prev_batch_index: 0,
            acc: r.acc.clone(),
            seq: r.seq.as_bytes().to_vec(),
            qual: r.qual.as_bytes().to_vec(),
            score: r.score,
        })
        .collect();

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
        let mut clusters = sweep::OrderedClusters::default();
        let mut reps: std::collections::HashMap<usize, sweep::ReadInfo> =
            std::collections::HashMap::new();
        for x in &sweep_reads {
            clusters.insert(x.id, vec![x.acc.clone()]);
            reps.insert(
                x.id,
                sweep::ReadInfo {
                    id: x.id,
                    batch_index: x.prev_batch_index,
                    acc: x.acc.clone(),
                    seq: x.seq.clone(),
                    qual: x.qual.clone(),
                    score: x.score,
                    error_rate: None,
                },
            );
        }
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

    let mut clusters_out = String::new();
    let mut origins_out = String::new();
    let mut nontrivial = 0usize;
    for (output_cl_id, c_id) in order.iter().enumerate() {
        let rep = &representatives[c_id];
        origins_out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            output_cl_id,
            strip_score(&rep.acc),
            String::from_utf8_lossy(&rep.seq),
            String::from_utf8_lossy(&rep.qual),
            pyfloat::repr(rep.score),
            pyfloat::repr(rep.error_rate.unwrap_or(f64::NAN)),
        ));
        let mut members: Vec<&String> = clusters.map[c_id].iter().collect();
        members.sort_by(|a, b| {
            score_of(b)
                .partial_cmp(&score_of(a))
                .expect("scores parse from the accession")
        });
        for acc in &members {
            clusters_out.push_str(&format!("{}\t{}\n", output_cl_id, strip_score(acc)));
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

    println!("Total number of reads iterated through:{}", reads.len());
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

    let ctext = match std::fs::read_to_string(clusters_path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("isONclust: cannot read {clusters_path}: {e}");
            return ExitCode::from(1);
        }
    };
    let mut order: Vec<String> = Vec::new();
    let mut members: std::collections::HashMap<String, Vec<String>> =
        std::collections::HashMap::new();
    for line in ctext.lines() {
        // `line.strip().split()` -- any whitespace, and blank lines vanish.
        let mut it = line.split_whitespace();
        let (cl_id, acc) = match (it.next(), it.next()) {
            (Some(a), Some(b)) => (a, b),
            _ => continue,
        };
        members
            .entry(cl_id.to_string())
            .or_insert_with(|| {
                order.push(cl_id.to_string());
                Vec::new()
            })
            .push(acc.to_string());
    }

    let ftext = match std::fs::read_to_string(fastq_path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("isONclust: cannot read {fastq_path}: {e}");
            return ExitCode::from(1);
        }
    };
    let mut reads: std::collections::HashMap<String, (String, String)> =
        std::collections::HashMap::new();
    for r in fastq::read(&ftext) {
        reads.insert(r.name, (r.seq, r.qual.unwrap_or_default()));
    }

    if let Err(e) = std::fs::create_dir_all(outfolder) {
        eprintln!("isONclust: cannot create {outfolder}: {e}");
        return ExitCode::from(1);
    }
    for cl_id in &order {
        let accs = &members[cl_id];
        if (accs.len() as i64) < wf.n {
            continue;
        }
        let mut body = String::new();
        for acc in accs {
            match reads.get(acc) {
                Some((seq, qual)) => {
                    body.push_str(&format!("@{acc}\n{seq}\n+\n{qual}\n"));
                }
                None => {
                    // The reference raises KeyError here.
                    eprintln!(
                        "isONclust: read {acc:?} is in {clusters_path} but not in {fastq_path}"
                    );
                    return ExitCode::from(1);
                }
            }
        }
        let path = std::path::Path::new(outfolder).join(format!("{cl_id}.fastq"));
        if let Err(e) = std::fs::write(&path, body) {
            eprintln!("isONclust: cannot write {}: {e}", path.display());
            return ExitCode::from(1);
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
