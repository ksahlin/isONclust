//! argparse-compatible command line parsing.
//!
//! This is hand-written rather than built on `clap`, and that is a deliberate
//! choice. The contract is the reference's *exact* stdout, stderr and exit code
//! (`bench/golden/cli/`, 28 cases), and argparse differs from clap on four
//! independent axes:
//!
//! 1. **Prefix abbreviation.** argparse accepts any unambiguous prefix, so
//!    `--outfold` works and `--f` is an error naming the candidates. clap has
//!    no equivalent.
//! 2. **Error text.** `isONclust: error: argument --k: invalid int value: 'abc'`
//!    is argparse's format, preceded by a fixed usage block, exit 2.
//! 3. **Double-dash single-letter options.** `--t`, `--d`, `--q`, `--k`, `--w`
//!    are long options one character wide, which is not what clap's `short`
//!    produces.
//! 4. **Exit codes.** Four validation paths exit **0** on what is logically an
//!    error (see `validate`).
//!
//! Bending clap to all four would have been more code than this, and every
//! bend is a place to diverge silently.

use crate::text;

/// What the program should do after parsing.
pub enum Outcome {
    /// Print to stdout and exit 0 (`--help`, `--version`, no arguments).
    Exit0(String),
    /// Print the usage block plus this line to stderr and exit 2.
    UsageError(String),
    /// Print to stdout and exit with this code (the validation paths).
    Message(String, i32),
    /// Parsed and validated; run the clustering.
    Run(Box<Args>),
    /// Parsed and validated; run the `write_fastq` subcommand.
    WriteFastq(WriteFastqArgs),
}

#[derive(Debug, Clone, PartialEq)]
pub struct Args {
    pub fastq: Option<String>,
    pub flnc: Option<String>,
    pub ccs: Option<String>,
    pub nr_cores: i64,
    pub print_output: i64,
    pub quality_threshold: f64,
    pub ont: bool,
    pub isoseq: bool,
    pub use_old_sorted_file: bool,
    pub consensus: bool,
    pub abundance_ratio: f64,
    pub rc_identity_threshold: f64,
    pub medaka: bool,
    pub k: i64,
    pub w: i64,
    pub min_shared: i64,
    pub mapped_threshold: f64,
    pub aligned_threshold: f64,
    pub batch_type: String,
    pub min_fraction: f64,
    pub min_prob_no_hits: f64,
    pub outfolder: Option<String>,
}

impl Default for Args {
    /// Exactly the reference's `argparse` defaults.
    fn default() -> Self {
        Args {
            fastq: None,
            flnc: None,
            ccs: None,
            nr_cores: 8,
            print_output: 10000,
            quality_threshold: 7.0,
            ont: false,
            isoseq: false,
            use_old_sorted_file: false,
            consensus: false,
            abundance_ratio: 0.1,
            rc_identity_threshold: 0.9,
            medaka: false,
            k: 15,
            w: 50,
            min_shared: 5,
            mapped_threshold: 0.7,
            aligned_threshold: 0.4,
            batch_type: "total_nt".to_string(),
            min_fraction: 0.8,
            min_prob_no_hits: 0.1,
            outfolder: None,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct WriteFastqArgs {
    pub clusters: Option<String>,
    pub fastq: Option<String>,
    pub outfolder: Option<String>,
    pub n: i64,
}

/// Option kinds, mirroring the `add_argument` calls one for one.
#[derive(Clone, Copy, PartialEq)]
enum Kind {
    Str,
    Int,
    Float,
    Flag,
}

struct Opt {
    name: &'static str,
    kind: Kind,
}

const fn o(name: &'static str, kind: Kind) -> Opt {
    Opt { name, kind }
}

/// The main parser's options, in `add_argument` order. `--help` and `--version`
/// are actions that fire during parsing and win over everything.
const MAIN_OPTS: &[Opt] = &[
    o("--help", Kind::Flag),
    o("--version", Kind::Flag),
    o("--fastq", Kind::Str),
    o("--flnc", Kind::Str),
    o("--ccs", Kind::Str),
    o("--t", Kind::Int),
    o("--d", Kind::Int),
    o("--q", Kind::Float),
    o("--ont", Kind::Flag),
    o("--isoseq", Kind::Flag),
    o("--use_old_sorted_file", Kind::Flag),
    o("--consensus", Kind::Flag),
    o("--abundance_ratio", Kind::Float),
    o("--rc_identity_threshold", Kind::Float),
    o("--medaka", Kind::Flag),
    o("--k", Kind::Int),
    o("--w", Kind::Int),
    o("--min_shared", Kind::Int),
    o("--mapped_threshold", Kind::Float),
    o("--aligned_threshold", Kind::Float),
    o("--batch_type", Kind::Str),
    o("--min_fraction", Kind::Float),
    o("--min_prob_no_hits", Kind::Float),
    o("--outfolder", Kind::Str),
];

const WF_OPTS: &[Opt] = &[
    o("--help", Kind::Flag),
    o("--clusters", Kind::Str),
    o("--fastq", Kind::Str),
    o("--outfolder", Kind::Str),
    o("--N", Kind::Int),
];

/// Resolve a token to a full option name, argparse-style.
///
/// An exact match always wins, even when the name is a prefix of others
/// (`--d` is an option and also a prefix of nothing, but `--k` and `--w` show
/// the shape). Otherwise a unique prefix is accepted, and an ambiguous one is
/// an error naming every candidate in the parser's declaration order --
/// argparse lists them in `add_argument` order, not sorted.
fn resolve<'a>(token: &str, opts: &'a [Opt]) -> Result<&'a Opt, String> {
    if let Some(x) = opts.iter().find(|x| x.name == token) {
        return Ok(x);
    }
    let matches: Vec<&Opt> = opts.iter().filter(|x| x.name.starts_with(token)).collect();
    match matches.len() {
        1 => Ok(matches[0]),
        0 => Err(String::new()), // caller reports "unrecognized arguments"
        _ => {
            let names: Vec<&str> = matches.iter().map(|x| x.name).collect();
            Err(format!(
                "ambiguous option: {} could match {}",
                token,
                names.join(", ")
            ))
        }
    }
}

fn type_error(name: &str, kind: Kind, raw: &str) -> String {
    let ty = match kind {
        Kind::Int => "int",
        Kind::Float => "float",
        _ => unreachable!("only Int and Float can fail to convert"),
    };
    format!("argument {}: invalid {} value: '{}'", name, ty, raw)
}

/// Python's `int()` on a CLI value: no underscores, no float syntax. `--t 1.5`
/// is an error, which is why this is not `str::parse::<i64>` alone -- it is,
/// but the point is that argparse does not accept `1.5` for an int and neither
/// does this.
fn parse_int(name: &str, raw: &str) -> Result<i64, String> {
    raw.trim()
        .parse::<i64>()
        .map_err(|_| type_error(name, Kind::Int, raw))
}

fn parse_float(name: &str, raw: &str) -> Result<f64, String> {
    raw.trim()
        .parse::<f64>()
        .map_err(|_| type_error(name, Kind::Float, raw))
}

/// Parse `argv` (without the program name).
pub fn parse(argv: &[String]) -> Outcome {
    // The `write_fastq` subcommand is a positional. argparse dispatches on the
    // first non-option token.
    if let Some(first) = argv.first() {
        if !first.starts_with('-') {
            if first == "write_fastq" {
                return parse_write_fastq(&argv[1..]);
            }
            return Outcome::UsageError(format!(
                "argument {{write_fastq}}: invalid choice: '{}' (choose from write_fastq)",
                first
            ));
        }
    }

    let mut args = Args::default();
    let mut unrecognized: Vec<String> = Vec::new();
    // Which of --k/--w the user set explicitly. Not used to resolve them (the
    // presets overwrite regardless -- see `validate`), but kept because it is
    // the sort of thing a later stage asks for and getting it wrong is silent.
    let mut i = 0usize;

    while i < argv.len() {
        let tok = &argv[i];

        if tok == "--" {
            i += 1;
            continue;
        }
        if !tok.starts_with('-') {
            unrecognized.push(tok.clone());
            i += 1;
            continue;
        }

        // `-h` is the only short option argparse defines here.
        if tok == "-h" {
            return Outcome::Exit0(text::HELP.to_string());
        }

        let (name_tok, inline) = match tok.split_once('=') {
            Some((n, v)) => (n.to_string(), Some(v.to_string())),
            None => (tok.clone(), None),
        };

        let opt = match resolve(&name_tok, MAIN_OPTS) {
            Ok(x) => x,
            Err(msg) if msg.is_empty() => {
                unrecognized.push(tok.clone());
                i += 1;
                continue;
            }
            Err(msg) => return Outcome::UsageError(msg),
        };

        // Actions fire immediately, during parsing, and win over everything --
        // including arguments that would otherwise be errors.
        match opt.name {
            "--help" => return Outcome::Exit0(text::HELP.to_string()),
            "--version" => return Outcome::Exit0(text::VERSION.to_string()),
            _ => {}
        }

        if opt.kind == Kind::Flag {
            match opt.name {
                "--ont" => args.ont = true,
                "--isoseq" => args.isoseq = true,
                "--use_old_sorted_file" => args.use_old_sorted_file = true,
                "--consensus" => args.consensus = true,
                "--medaka" => args.medaka = true,
                _ => unreachable!("unhandled flag {}", opt.name),
            }
            i += 1;
            continue;
        }

        let raw = match inline {
            Some(v) => v,
            None => {
                i += 1;
                match argv.get(i) {
                    Some(v) => v.clone(),
                    None => {
                        return Outcome::UsageError(format!(
                            "argument {}: expected one argument",
                            opt.name
                        ))
                    }
                }
            }
        };

        let set = |args: &mut Args| -> Result<(), String> {
            match opt.name {
                "--fastq" => args.fastq = Some(raw.clone()),
                "--flnc" => args.flnc = Some(raw.clone()),
                "--ccs" => args.ccs = Some(raw.clone()),
                "--batch_type" => args.batch_type = raw.clone(),
                "--outfolder" => args.outfolder = Some(raw.clone()),
                "--t" => args.nr_cores = parse_int(opt.name, &raw)?,
                "--d" => args.print_output = parse_int(opt.name, &raw)?,
                "--k" => args.k = parse_int(opt.name, &raw)?,
                "--w" => args.w = parse_int(opt.name, &raw)?,
                "--min_shared" => args.min_shared = parse_int(opt.name, &raw)?,
                "--q" => args.quality_threshold = parse_float(opt.name, &raw)?,
                "--abundance_ratio" => args.abundance_ratio = parse_float(opt.name, &raw)?,
                "--rc_identity_threshold" => {
                    args.rc_identity_threshold = parse_float(opt.name, &raw)?
                }
                "--mapped_threshold" => args.mapped_threshold = parse_float(opt.name, &raw)?,
                "--aligned_threshold" => args.aligned_threshold = parse_float(opt.name, &raw)?,
                "--min_fraction" => args.min_fraction = parse_float(opt.name, &raw)?,
                "--min_prob_no_hits" => args.min_prob_no_hits = parse_float(opt.name, &raw)?,
                _ => unreachable!("unhandled option {}", opt.name),
            }
            Ok(())
        };
        if let Err(msg) = set(&mut args) {
            return Outcome::UsageError(msg);
        }
        i += 1;
    }

    if !unrecognized.is_empty() {
        return Outcome::UsageError(format!(
            "unrecognized arguments: {}",
            unrecognized.join(" ")
        ));
    }

    validate(args, argv.is_empty())
}

fn parse_write_fastq(argv: &[String]) -> Outcome {
    let mut wf = WriteFastqArgs::default();
    let mut unrecognized: Vec<String> = Vec::new();
    let mut i = 0usize;
    while i < argv.len() {
        let tok = &argv[i];
        if tok == "-h" {
            return Outcome::Exit0(text::WRITE_FASTQ_HELP.to_string());
        }
        if !tok.starts_with('-') {
            unrecognized.push(tok.clone());
            i += 1;
            continue;
        }
        let (name_tok, inline) = match tok.split_once('=') {
            Some((n, v)) => (n.to_string(), Some(v.to_string())),
            None => (tok.clone(), None),
        };
        let opt = match resolve(&name_tok, WF_OPTS) {
            Ok(x) => x,
            Err(msg) if msg.is_empty() => {
                unrecognized.push(tok.clone());
                i += 1;
                continue;
            }
            Err(msg) => return Outcome::UsageError(msg),
        };
        if opt.name == "--help" {
            return Outcome::Exit0(text::WRITE_FASTQ_HELP.to_string());
        }
        let raw = match inline {
            Some(v) => v,
            None => {
                i += 1;
                match argv.get(i) {
                    Some(v) => v.clone(),
                    None => {
                        return Outcome::UsageError(format!(
                            "argument {}: expected one argument",
                            opt.name
                        ))
                    }
                }
            }
        };
        match opt.name {
            "--clusters" => wf.clusters = Some(raw),
            "--fastq" => wf.fastq = Some(raw),
            "--outfolder" => wf.outfolder = Some(raw),
            "--N" => match parse_int(opt.name, &raw) {
                Ok(v) => wf.n = v,
                Err(msg) => return Outcome::UsageError(msg),
            },
            _ => unreachable!("unhandled write_fastq option {}", opt.name),
        }
        i += 1;
    }
    if !unrecognized.is_empty() {
        // argparse reports stray tokens against the MAIN parser's usage, because
        // the subparser hands them back up. Measured: `write_fastq --k 5` gives
        // "unrecognized arguments: --k 5".
        return Outcome::UsageError(format!(
            "unrecognized arguments: {}",
            unrecognized.join(" ")
        ));
    }
    Outcome::WriteFastq(wf)
}

/// Post-parse validation, in the reference's exact order.
///
/// Four of these exit **0** on what is logically an error, because the
/// reference calls bare `sys.exit()`. That is the contract, warts included: a
/// wrapper script cannot currently detect them. Changing it is in
/// PORTING.md's *Deferred improvements*, not here.
fn validate(mut args: Args, argv_was_empty: bool) -> Outcome {
    let has = |o: &Option<String>| o.is_some();

    // 1. fastq together with either BAM input
    if has(&args.fastq) && (has(&args.flnc) || has(&args.ccs)) {
        return Outcome::Message(
            "Either (1) only a fastq file, or (2) a ccs and a flnc file should be specified. \n"
                .to_string(),
            0,
        );
    }
    // 2. exactly one of flnc/ccs
    if has(&args.flnc) != has(&args.ccs) {
        return Outcome::Message(
            "isONclust needs both the ccs.bam file produced by ccs and the flnc file produced by isoseq3 cluster. \n"
                .to_string(),
            0,
        );
    }
    // 3. the presets. Note these OVERWRITE an explicit --k/--w regardless of
    //    the order they appeared in: `--ont --k 99` resolves to k=13, and so
    //    does `--k 99 --ont`. Measured; see bench/golden/cli/ont_over_k.
    if args.ont && args.isoseq {
        return Outcome::Message(
            "Arguments mutually exclusive, specify either --isoseq or --ont. \n".to_string(),
            0,
        );
    } else if args.isoseq {
        args.k = 15;
        args.w = 50;
    } else if args.ont {
        args.k = 13;
        args.w = 20;
    }
    // 4. no arguments at all
    if argv_was_empty {
        return Outcome::Exit0(text::HELP.to_string());
    }
    // 5. no input specified
    if !has(&args.fastq) && !has(&args.flnc) && !has(&args.ccs) {
        return Outcome::Exit0(text::HELP.to_string());
    }
    // 6. window size. The only validation that exits non-zero. One message for
    //    both directions, and note `100 < w` is checked before `w < k`.
    if 100 < args.w || args.w < args.k {
        return Outcome::Message(
            "Please specify a window of size larger or equal to k, and smaller than 100.\n"
                .to_string(),
            1,
        );
    }
    Outcome::Run(Box::new(args))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(args: &[&str]) -> Vec<String> {
        args.iter().map(|s| s.to_string()).collect()
    }

    /// Parse and expect a successful run, returning the resolved args.
    fn run(args: &[&str]) -> Args {
        match parse(&v(args)) {
            Outcome::Run(a) => *a,
            _ => panic!("expected Run for {:?}", args),
        }
    }

    fn usage_err(args: &[&str]) -> String {
        match parse(&v(args)) {
            Outcome::UsageError(m) => m,
            _ => panic!("expected UsageError for {:?}", args),
        }
    }

    fn message(args: &[&str]) -> (String, i32) {
        match parse(&v(args)) {
            Outcome::Message(m, c) => (m, c),
            _ => panic!("expected Message for {:?}", args),
        }
    }

    #[test]
    fn defaults_match_argparse() {
        let a = run(&["--fastq", "x"]);
        assert_eq!(a.k, 15);
        assert_eq!(a.w, 50);
        assert_eq!(a.nr_cores, 8);
        assert_eq!(a.print_output, 10000);
        assert_eq!(a.quality_threshold, 7.0);
        assert_eq!(a.min_shared, 5);
        assert_eq!(a.mapped_threshold, 0.7);
        assert_eq!(a.aligned_threshold, 0.4);
        assert_eq!(a.min_fraction, 0.8);
        assert_eq!(a.min_prob_no_hits, 0.1);
        assert_eq!(a.batch_type, "total_nt");
        assert_eq!(a.abundance_ratio, 0.1);
        assert_eq!(a.rc_identity_threshold, 0.9);
        assert!(!a.ont && !a.isoseq && !a.consensus && !a.medaka);
    }

    #[test]
    fn presets_are_exactly_k_and_w() {
        let a = run(&["--ont", "--fastq", "x"]);
        assert_eq!((a.k, a.w), (13, 20));
        let b = run(&["--isoseq", "--fastq", "x"]);
        assert_eq!((b.k, b.w), (15, 50));
    }

    /// The wart worth pinning: a preset silently discards an explicit --k/--w,
    /// whichever order they appear in. Measured against the reference; see
    /// bench/golden/cli/ont_over_k and k_then_ont.
    #[test]
    fn preset_overwrites_explicit_k_and_w_in_either_order() {
        assert_eq!(run(&["--ont", "--k", "99", "--fastq", "x"]).k, 13);
        assert_eq!(run(&["--k", "99", "--ont", "--fastq", "x"]).k, 13);
        assert_eq!(run(&["--isoseq", "--w", "7", "--fastq", "x"]).w, 50);
        assert_eq!(run(&["--w", "7", "--isoseq", "--fastq", "x"]).w, 50);
    }

    #[test]
    fn unique_prefix_is_accepted() {
        assert_eq!(run(&["--fas", "x"]).fastq.as_deref(), Some("x"));
        assert_eq!(
            run(&["--fastq", "x", "--outfold", "o"])
                .outfolder
                .as_deref(),
            Some("o")
        );
        assert_eq!(run(&["--fastq", "x", "--min_sh", "9"]).min_shared, 9);
    }

    #[test]
    fn ambiguous_prefix_names_every_candidate_in_declaration_order() {
        assert_eq!(
            usage_err(&["--f", "x"]),
            "ambiguous option: --f could match --fastq, --flnc"
        );
        assert_eq!(
            usage_err(&["--m", "5", "--fastq", "x"]),
            "ambiguous option: --m could match --medaka, --min_shared, --mapped_threshold, --min_fraction, --min_prob_no_hits"
        );
    }

    /// An exact name wins even when it is a prefix of longer options.
    #[test]
    fn exact_match_beats_prefix() {
        assert_eq!(run(&["--fastq", "x", "--k", "9", "--w", "20"]).k, 9);
        assert_eq!(run(&["--fastq", "x", "--d", "5"]).print_output, 5);
    }

    #[test]
    fn type_errors_use_argparse_wording() {
        assert_eq!(
            usage_err(&["--k", "abc", "--fastq", "x"]),
            "argument --k: invalid int value: 'abc'"
        );
        assert_eq!(
            usage_err(&["--q", "xyz", "--fastq", "x"]),
            "argument --q: invalid float value: 'xyz'"
        );
        // an int option rejects float syntax
        assert_eq!(
            usage_err(&["--t", "1.5", "--fastq", "x"]),
            "argument --t: invalid int value: '1.5'"
        );
    }

    #[test]
    fn missing_value_is_reported_against_the_full_option_name() {
        assert_eq!(
            usage_err(&["--fastq", "x", "--k"]),
            "argument --k: expected one argument"
        );
        // ...even when the user typed a prefix
        assert_eq!(
            usage_err(&["--fastq", "x", "--min_sh"]),
            "argument --min_shared: expected one argument"
        );
    }

    #[test]
    fn equals_form_is_accepted() {
        assert_eq!(run(&["--fastq=x", "--k=9", "--w=20"]).k, 9);
        assert_eq!(run(&["--fastq=x"]).fastq.as_deref(), Some("x"));
    }

    #[test]
    fn unrecognized_arguments_are_collected_and_joined() {
        assert_eq!(
            usage_err(&["--fastq", "x", "--no-such-flag"]),
            "unrecognized arguments: --no-such-flag"
        );
    }

    #[test]
    fn help_and_version_win_over_broken_arguments() {
        assert!(matches!(
            parse(&v(&["--version", "--fastq", "/nope", "--k", "abc"])),
            Outcome::Exit0(_)
        ));
        assert!(matches!(
            parse(&v(&["--fastq", "/nope", "-h"])),
            Outcome::Exit0(_)
        ));
        match parse(&v(&["--version"])) {
            Outcome::Exit0(s) => assert_eq!(s, "isONclust 0.0.6.1\n"),
            _ => panic!("expected Exit0"),
        }
    }

    /// Four validation paths exit 0 on what is logically an error. This is the
    /// reference's behaviour and it is deliberately reproduced; see the note on
    /// `validate`.
    #[test]
    fn validation_paths_that_exit_zero() {
        for args in [
            vec!["--fastq", "x", "--ccs", "y"],
            vec!["--flnc", "y"],
            vec!["--ccs", "y"],
            vec!["--ont", "--isoseq", "--fastq", "x"],
        ] {
            let (_, code) = message(&args);
            assert_eq!(code, 0, "{:?} should exit 0", args);
        }
    }

    #[test]
    fn window_validation_is_the_only_nonzero_exit() {
        let (msg, code) = message(&["--fastq", "x", "--k", "20", "--w", "15"]);
        assert_eq!(code, 1);
        assert!(msg.starts_with("Please specify a window"));
        // both directions share one message
        let (msg2, code2) = message(&["--fastq", "x", "--k", "15", "--w", "101"]);
        assert_eq!(code2, 1);
        assert_eq!(msg, msg2);
        // w == k is allowed
        assert!(matches!(
            parse(&v(&["--fastq", "x", "--k", "15", "--w", "15"])),
            Outcome::Run(_)
        ));
    }

    #[test]
    fn no_arguments_and_no_input_both_print_help() {
        assert!(matches!(parse(&[]), Outcome::Exit0(_)));
        assert!(matches!(parse(&v(&["--k", "9"])), Outcome::Exit0(_)));
    }

    #[test]
    fn write_fastq_subcommand() {
        match parse(&v(&[
            "write_fastq",
            "--clusters",
            "c",
            "--fastq",
            "f",
            "--outfolder",
            "o",
            "--N",
            "2",
        ])) {
            Outcome::WriteFastq(w) => {
                assert_eq!(w.clusters.as_deref(), Some("c"));
                assert_eq!(w.fastq.as_deref(), Some("f"));
                assert_eq!(w.outfolder.as_deref(), Some("o"));
                assert_eq!(w.n, 2);
            }
            _ => panic!("expected WriteFastq"),
        }
        // default --N
        match parse(&v(&["write_fastq"])) {
            Outcome::WriteFastq(w) => assert_eq!(w.n, 0),
            _ => panic!("expected WriteFastq"),
        }
    }

    #[test]
    fn write_fastq_rejects_main_parser_flags() {
        assert_eq!(
            usage_err(&["write_fastq", "--k", "5"]),
            "unrecognized arguments: --k 5"
        );
    }

    #[test]
    fn unknown_subcommand_lists_the_choices() {
        assert_eq!(
            usage_err(&["bogus"]),
            "argument {write_fastq}: invalid choice: 'bogus' (choose from write_fastq)"
        );
    }
}
