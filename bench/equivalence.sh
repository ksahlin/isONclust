#!/usr/bin/env bash
#
# The equivalence harness. Equivalence with the Python reference is the
# acceptance criterion for the port, and this is what checks it.
#
#   bench/equivalence.sh env      # is the reference environment usable?
#   bench/equivalence.sh seeds    # determinism gate -- run this BEFORE recording
#   bench/equivalence.sh cli      # CLI contract: exit codes, stdout, stderr
#   bench/equivalence.sh record   # record goldens from the reference
#   bench/equivalence.sh verify   # run the port, diff against the goldens
#   bench/equivalence.sh dropped  # the port must refuse the out-of-scope flags
#   bench/equivalence.sh stable   # recording twice must give identical goldens
#   bench/equivalence.sh all      # everything
#
# Environment:
#   REF_PYTHON   interpreter that can import parasail and pysam
#                (default: the isonclust-ref conda env; see setup_reference_env.sh)
#   PORT_BIN     the Rust binary under test (default: rust/target/release/isONclust)
#   CORPUS       input fastq path, or a name from bench/corpora.tsv (default: smoke)
#   GOLDEN       where goldens live (default: bench/golden)
#
# WHAT COUNTS AS A DIFFERENCE
# ---------------------------
# Every file the tool writes, byte for byte:
#   final_clusters.tsv          the actual result
#   final_cluster_origins.tsv   representative per cluster, and its error_rate
#   sorted.fastq                the sorted input, with the score in each header
#   logfile.txt                 error-rate summary statistics
#   <n>/pre_clusters.csv        parallel mode only, one dir per merge iteration
#   <n>/cluster_origins.csv     parallel mode only
#
# sorted.fastq is NOT an intermediate to be skipped. The score is formatted into
# every read accession with Python's float repr and then parsed back out with
# float(), so it is both an output and an input, and Rust's default float
# formatting does not match Python's (Python writes 1234.0 and 1e-05 where Rust
# writes 1234 and 0.00001).
#
# logfile.txt is included because it is the only place the error-rate
# distribution is observable, and error_rate is where the reference's one
# determinism defect shows up. See `seeds`.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
cd "$ROOT"

REF_PYTHON="${REF_PYTHON:-$HOME/miniforge3/envs/isonclust-ref/bin/python}"
PORT_BIN="${PORT_BIN:-$ROOT/rust/target/release/isONclust}"
# CORPUS accepts a path, or a name from bench/corpora.tsv.
ISONCLUST_DATA="${ISONCLUST_DATA:-$HOME/data/lrRNA-seq}"
resolve_corpus() {
  local c="$1" p
  [[ -f "$c" ]] && { echo "$c"; return; }
  p="$(awk -F'\t' -v n="$c" '$1==n {print $2; exit}' "$ROOT/bench/corpora.tsv" 2>/dev/null)"
  [[ -z "$p" ]] && { echo "$c"; return; }          # not a known name: pass through and let it fail visibly
  case "$p" in
    /*)     echo "$p" ;;
    test/*) echo "$ROOT/$p" ;;
    *)      echo "$ISONCLUST_DATA/$p" ;;
  esac
}
CORPUS="$(resolve_corpus "${CORPUS:-smoke}")"
GOLDEN="${GOLDEN:-$ROOT/bench/golden}"
WORK="${WORK:-$(mktemp -d)}"

OUT_FILES=(final_clusters.tsv final_cluster_origins.tsv sorted.fastq logfile.txt)

PASS=0; FAIL=0
ok()   { printf '    ok    %s\n' "$1"; PASS=$((PASS+1)); }
bad()  { printf '    FAIL  %s\n' "$1"; FAIL=$((FAIL+1)); }
info() { printf '    info  %s\n' "$1"; }

# ---------------------------------------------------------------------------

cmd_env() {
  echo "==> reference environment"
  [[ -x "$REF_PYTHON" ]] || { bad "REF_PYTHON not executable: $REF_PYTHON"; return; }
  if "$REF_PYTHON" - <<'PY'
import sys
import parasail, pysam
m = parasail.matrix_create("ACGT", 2, -2)
r = parasail.sg_trace_scan_16("ACGTACGTAA", "ACGTTCGTAA", 5, 1, m)
assert r.score == 16 and str(r.cigar.decode, "utf-8") == "4=1X5="
xs = [0.1] * 10 + [1e17, -1e17]
compensated = sum(xs) == sum(reversed(xs))
print(f"    ok    python {sys.version.split()[0]}, pysam {pysam.__version__}, parasail 1.3.4")
print(f"    ok    sum() is {'COMPENSATED (>=3.12)' if compensated else 'NAIVE (<=3.11)'}")
if not compensated:
    print("    info  on this interpreter the reference is NOT deterministic run to run;")
    print("    info  `equivalence.sh seeds` will fail. See PORTING.md.")
PY
  then PASS=$((PASS+2)); else bad "reference env unusable"; fi
  [[ -f "$CORPUS" ]] && ok "corpus present: $CORPUS ($(( $(wc -l < "$CORPUS") / 4 )) reads)" \
                     || bad "corpus missing: $CORPUS"
  # spoa is deliberately NOT checked. It was needed only by --consensus, which is
  # out of scope (PORTING.md, Scope), so the port needs no POA at all. The
  # reference environment is parasail + pysam and nothing else.
}

# ---------------------------------------------------------------------------
# The determinism gate. Recording a golden from a non-deterministic reference
# records one of its several possible answers, so this runs first and refuses to
# be skipped.

cmd_seeds() {
  echo "==> determinism: does the reference agree with itself across PYTHONHASHSEED?"
  local seeds=(0 1 2 7 12345) s d first
  for entry_args in "--isoseq --t 1" "--isoseq --t 8"; do
    first=""
    for s in "${seeds[@]}"; do
      d="$WORK/seeds_${s}_$(echo "$entry_args" | tr -dc 'a-z0-9')"
      rm -rf "$d"; mkdir -p "$d"
      PYTHONHASHSEED="$s" "$REF_PYTHON" isONclust $entry_args \
        --fastq "$CORPUS" --outfolder "$d" >"$d.stdout" 2>"$d.stderr" || true
    done
    for f in "${OUT_FILES[@]}"; do
      first=""; local differs=0
      for s in "${seeds[@]}"; do
        d="$WORK/seeds_${s}_$(echo "$entry_args" | tr -dc 'a-z0-9')"
        [[ -f "$d/$f" ]] || continue
        local h; h="$(shasum -a 256 "$d/$f" | cut -d' ' -f1)"
        [[ -z "$first" ]] && first="$h" && continue
        [[ "$h" == "$first" ]] || differs=1
      done
      if [[ $differs -eq 0 ]]; then
        ok "stable across ${#seeds[@]} seeds: $f  [$entry_args]"
      else
        bad "SEED-DEPENDENT: $f  [$entry_args]"
        # Say where, not just that. This is the actionable half.
        local a b
        a="$WORK/seeds_0_$(echo "$entry_args" | tr -dc 'a-z0-9')/$f"
        b="$WORK/seeds_1_$(echo "$entry_args" | tr -dc 'a-z0-9')/$f"
        info "  seed 0 vs seed 1:"
        "$REF_PYTHON" bench/diffsummary.py "$a" "$b" "$f" 2>/dev/null || true
      fi
    done
  done
}

# ---------------------------------------------------------------------------
# The CLI contract. Argument names, defaults, validation order, the exact stderr
# and stdout text, and the exit code. Recorded from the reference, then replayed
# against the port.
#
# Several of these exit 0 on what is logically an error -- no arguments, and
# `--ont --isoseq` together, both print a message and `sys.exit()` with no
# argument. That is the contract, warts included; do not "fix" it in the port
# without giving the divergence its own commit and a note in PORTING.md.

cli_case() { # cli_case <name> <args...>
  local name="$1"; shift
  local d="$GOLDEN/cli/$name"
  mkdir -p "$d"
  if [[ "$MODE" == "record" ]]; then
    set +e
    "$REF_PYTHON" isONclust "$@" >"$d/stdout" 2>"$d/stderr"; echo $? >"$d/exit"
    set -e
    # Paths and timings are not contract; scrub them so the golden is portable.
    sed -i.bak -E 's#(/private)?(/var/folders/[^ ]*|/tmp/[^ ]*)#<TMPDIR>#g; s#/[^ ]*/(isONclust|sirv_sim_120)#<PATH>/\1#g; s/[0-9]+\.[0-9]{4,}/<TIME>/g' "$d/stdout" "$d/stderr"
    rm -f "$d"/*.bak
    ok "recorded cli/$name (exit $(cat "$d/exit"))"
  else
    [[ -f "$d/exit" ]] || { bad "no golden for cli/$name -- re-run: equivalence.sh cli record"; return; }
    set +e
    "$PORT_BIN" "$@" >"$WORK/o" 2>"$WORK/e"; local rc=$?
    set -e
    sed -i.bak -E 's#(/private)?(/var/folders/[^ ]*|/tmp/[^ ]*)#<TMPDIR>#g; s#/[^ ]*/(isONclust|sirv_sim_120)#<PATH>/\1#g; s/[0-9]+\.[0-9]{4,}/<TIME>/g' "$WORK/o" "$WORK/e"
    if [[ "$rc" == "$(cat "$d/exit")" ]] && diff -q "$WORK/e" "$d/stderr" >/dev/null; then
      ok "cli/$name"
    else
      bad "cli/$name (exit $rc, want $(cat "$d/exit"))"
      diff "$d/stderr" "$WORK/e" | head -6 | sed 's/^/          /'
    fi
  fi
}

cmd_cli() {
  MODE="${1:-record}"
  echo "==> CLI contract ($MODE)"
  # Without this, verify mode returns from every case silently and the run
  # reports "0 passed, 0 failed", which reads like success.
  if [[ "$MODE" == "verify" && ! -x "$PORT_BIN" ]]; then
    bad "no port binary at $PORT_BIN -- nothing to verify yet"
    return
  fi
  cli_case version      --version
  cli_case help         --help
  cli_case noargs
  cli_case both_presets --ont --isoseq --fastq "$CORPUS"
  cli_case w_lt_k       --fastq "$CORPUS" --k 20 --w 15
  cli_case w_gt_100     --fastq "$CORPUS" --k 15 --w 101
  cli_case flnc_no_ccs  --flnc x.bam
  cli_case ccs_no_flnc  --ccs x.bam
  cli_case fastq_and_ccs --fastq "$CORPUS" --ccs x.bam
  cli_case unknown_flag --fastq "$CORPUS" --no-such-flag
  cli_case d_zero       --fastq "$CORPUS" --outfolder "$WORK/dz" --t 1 --d 0
  cli_case medaka       --fastq "$CORPUS" --outfolder "$WORK/mk" --t 1 --consensus --medaka
  # argparse accepts any unambiguous prefix; clap does not, and every
  # multi-word flag needs an explicit long name or clap renames it.
  cli_case abbrev_outf  --fastq "$CORPUS" --outfold "$WORK/ab" --t 1 --k 99
  cli_case wf_help      write_fastq --help
}

# ---------------------------------------------------------------------------

run_case() { # run_case <name> <entry> <args>  -> populates $2 output dir
  local name="$1" entry="$2" args="$3" runner="$4" outdir="$5"
  rm -rf "$outdir"; mkdir -p "$outdir"
  set +e
  if [[ "$entry" == "write_fastq" ]]; then
    # write_fastq consumes a clustering, so it needs one to exist first.
    local clusters="$GOLDEN/sample/final_clusters.tsv"
    [[ -f "$clusters" ]] || { echo "SKIP"; return; }
    $runner write_fastq --clusters "$clusters" --fastq "$CORPUS" \
            --outfolder "$outdir" $args >"$outdir.stdout" 2>"$outdir.stderr"
  else
    $runner $args --fastq "$CORPUS" --outfolder "$outdir" \
            >"$outdir.stdout" 2>"$outdir.stderr"
  fi
  local rc=$?
  set -e
  echo "$rc"
}


# A case line with no tab silently degrades into "run with no arguments", which
# then records a plausible-looking golden for the wrong invocation. Refuse.
check_cases() {
  local bad_lines
  bad_lines="$(grep -vE '^#|^$' bench/cases.tsv | grep -vcE $'^[^\t]+\t[^\t]+\t' || true)"
  if [[ "$bad_lines" != "0" ]]; then
    bad "bench/cases.tsv has $bad_lines line(s) that are not TAB-separated into 3 fields"
    info "an editor probably expanded tabs to spaces; see the comment in that file"
    exit 1
  fi
  info "$(grep -vcE '^#|^$' bench/cases.tsv) cases, all tab-separated"
}

cmd_record() {
  echo "==> recording goldens from the reference"
  check_cases
  mkdir -p "$GOLDEN/out"

  # Goldens are a MANIFEST OF HASHES, not the files themselves. Recorded
  # verbatim, the 27 cases come to 318 MB: final_cluster_origins.tsv carries
  # every representative's full sequence and quality string (3.5 MB a case),
  # sorted.fastq is the whole input again (6.5 MB a case), and `wf_N0` alone is
  # 1240 files. None of that belongs in a repository this exercise just shrank
  # from 492 MB to 1 MB.
  #
  # Hashes are enough to FAIL correctly. They are not enough to say what broke,
  # so `verify` re-runs the reference for a failing case and diffs against that.
  # The reference is 2 seconds a case; the storage is not worth it.
  {
    echo "# isONclust reference goldens -- per-file sha256"
    echo "# No timestamp on purpose: it made this file differ on every re-record,"
    echo "# which leaves git permanently dirty and hides real changes in the noise."
    echo "# git records when it was committed; what matters for validity is below."
    echo "# corpus:   $(basename "$CORPUS")  sha256 $(shasum -a 256 "$CORPUS" | cut -d' ' -f1)"
    "$REF_PYTHON" -c "import sys,parasail,pysam; print(f'# reference: python {sys.version.split()[0]}, pysam {pysam.__version__}, parasail 1.3.4')"
    "$REF_PYTHON" -c "xs=[0.1]*10+[1e17,-1e17]; print('# sum():    ' + ('compensated (>=3.12)' if sum(xs)==sum(reversed(xs)) else 'NAIVE (<=3.11) -- these goldens are NOT reproducible'))"
    echo "# PYTHONHASHSEED=0 for every case"
    echo "#"
    echo "# case	exit	relpath	sha256	bytes"
  } > "$GOLDEN/manifest.tsv"

  local n=0
  while IFS=$'\t' read -r name entry args; do
    [[ "$name" =~ ^# ]] && continue
    [[ -z "${name// }" ]] && continue
    local d="$WORK/rec/$name"
    local rc; rc="$(PYTHONHASHSEED=0 run_case "$name" "$entry" "$args" "$REF_PYTHON isONclust" "$d")"
    if [[ "$rc" == "SKIP" ]]; then info "skipped $name (needs the default case first)"; continue; fi
    local nf=0
    while IFS= read -r rel; do
      printf '%s\t%s\t%s\t%s\t%s\n' "$name" "$rc" "$rel" \
        "$(shasum -a 256 "$d/$rel" | cut -d' ' -f1)" "$(wc -c < "$d/$rel" | tr -d ' ')" \
        >> "$GOLDEN/manifest.tsv"
      nf=$((nf+1))
    done < <(cd "$d" && find . -type f | sed 's|^\./||' | LC_ALL=C sort)
    # Keep ONE case's small files verbatim, so there is something to read by eye
    # without running anything. sorted.fastq and final_cluster_origins.tsv are
    # excluded by size; their hashes are in the manifest like everything else.
    if [[ "$name" == "default" ]]; then
      mkdir -p "$GOLDEN/sample"
      cp "$d/final_clusters.tsv" "$d/logfile.txt" "$GOLDEN/sample/" 2>/dev/null || true
      head -8 "$d/sorted.fastq" > "$GOLDEN/sample/sorted.fastq.head" 2>/dev/null || true
      cut -f1,5,6 "$d/final_cluster_origins.tsv" > "$GOLDEN/sample/final_cluster_origins.id_score_errorrate.tsv" 2>/dev/null || true
    fi
    n=$((n+1))
    ok "recorded $name (exit $rc, $nf files)"
  done < bench/cases.tsv
  info "$n cases -> $GOLDEN/manifest.tsv ($(wc -c < "$GOLDEN/manifest.tsv" | tr -d ' ') bytes)"
}

cmd_verify() {
  echo "==> verifying the port against the goldens"
  check_cases
  [[ -f "$GOLDEN/manifest.tsv" ]] || { bad "no manifest -- run: bench/equivalence.sh record"; return; }
  if [[ ! -x "$PORT_BIN" ]]; then
    bad "no port binary at $PORT_BIN -- nothing to verify yet"
    info "expected until rust/ exists; build with"
    info "  cargo build --release --manifest-path rust/Cargo.toml"
    return
  fi
  while IFS=$'\t' read -r name entry args; do
    [[ "$name" =~ ^# ]] && continue
    [[ -z "${name// }" ]] && continue
    local d="$WORK/port/$name"
    local want_exit; want_exit="$(awk -F'\t' -v n="$name" '$1==n {print $2; exit}' "$GOLDEN/manifest.tsv")"
    [[ -n "$want_exit" ]] || { info "no golden for $name"; continue; }
    local rc; rc="$(run_case "$name" "$entry" "$args" "$PORT_BIN" "$d")"

    local mismatched=() missing=()
    while IFS=$'\t' read -r rel want_sha want_bytes; do
      if [[ ! -f "$d/$rel" ]]; then missing+=("$rel"); continue; fi
      local got; got="$(shasum -a 256 "$d/$rel" | cut -d' ' -f1)"
      [[ "$got" == "$want_sha" ]] || mismatched+=("$rel")
    done < <(awk -F'\t' -v n="$name" '$1==n {print $3"\t"$4"\t"$5}' "$GOLDEN/manifest.tsv")

    # A file the port writes that the reference does not is also a failure.
    local extra=()
    while IFS= read -r rel; do
      awk -F'\t' -v n="$name" -v r="$rel" '$1==n && $3==r {found=1} END {exit !found}' \
        "$GOLDEN/manifest.tsv" || extra+=("$rel")
    done < <(cd "$d" 2>/dev/null && find . -type f | sed 's|^\./||' | LC_ALL=C sort)

    if [[ ${#mismatched[@]} -eq 0 && ${#missing[@]} -eq 0 && ${#extra[@]} -eq 0 && "$rc" == "$want_exit" ]]; then
      ok "$name"
      continue
    fi
    bad "$name (exit $rc, want $want_exit)"
    [[ ${#missing[@]}  -gt 0 ]] && info "  not written by the port: ${missing[*]}"
    [[ ${#extra[@]}    -gt 0 ]] && info "  written by the port only: ${extra[*]}"
    # Hashes cannot say what moved, so re-derive the reference output for this
    # one case and diff properly. Two seconds beats 318 MB in git.
    if [[ ${#mismatched[@]} -gt 0 ]]; then
      info "  differing: ${mismatched[*]}"
      local r="$WORK/refre/$name"
      PYTHONHASHSEED=0 run_case "$name" "$entry" "$args" "$REF_PYTHON isONclust" "$r" >/dev/null
      for rel in "${mismatched[@]}"; do
        [[ -f "$r/$rel" && -f "$d/$rel" ]] || continue
        info "  --- $rel ---"
        "$REF_PYTHON" bench/diffsummary.py "$r/$rel" "$d/$rel" "$(basename "$rel")" || true
      done
    fi
  done < bench/cases.tsv
}


# ---------------------------------------------------------------------------
# Dropped flags. --consensus, --abundance_ratio, --rc_identity_threshold and
# --medaka are out of scope (PORTING.md, Scope) and the port must refuse them
# rather than accept-and-ignore: a pipeline that passes --consensus, gets exit 0
# and no consensus_references.fasta is worse off than one that fails.
#
# This CANNOT be a recorded golden, because the reference does not do it -- given
# spoa the reference happily runs --consensus. It is the port's own contract, so
# it is asserted directly. Two things must hold: non-zero exit, and the flag
# named somewhere in the output (a bare "unrecognised argument" is not enough to
# act on).
DROPPED_FLAGS=(--consensus --medaka --abundance_ratio --rc_identity_threshold)

cmd_dropped() {
  echo "==> dropped flags must be refused, by name"
  if [[ ! -x "$PORT_BIN" ]]; then
    bad "no port binary at $PORT_BIN -- nothing to check yet"
    return
  fi
  local f arg rc out
  for f in "${DROPPED_FLAGS[@]}"; do
    case "$f" in
      --abundance_ratio|--rc_identity_threshold) arg="0.5" ;;   # these take a value
      *) arg="" ;;
    esac
    set +e
    out="$("$PORT_BIN" --fastq "$CORPUS" --outfolder "$WORK/dropped" --t 1 $f $arg 2>&1)"
    rc=$?
    set -e
    local name="${f#--}"
    if [[ $rc -eq 0 ]]; then
      bad "$f accepted (exit 0) -- must be refused"
    elif ! grep -qF -- "$name" <<<"$out"; then
      bad "$f refused (exit $rc) but the message does not name it"
      sed 's/^/          /' <<<"$out" | head -3
    else
      ok "$f refused, exit $rc, named in the message"
    fi
  done
}

# ---------------------------------------------------------------------------
# Recording twice must give the same goldens. A golden containing a timestamp, a
# temp path, a PID or a duration can never be matched by anything -- including
# the reference itself -- so it is not a check, it is a permanent failure that
# trains you to ignore the harness. Two got through: a `tempfile.mkdtemp()` path
# in the --consensus stdout, and the manifest's own "recorded:" line.

cmd_stable() {
  echo "==> recording twice must be byte-identical"
  local a="$WORK/stable_a" b="$WORK/stable_b"
  GOLDEN="$a" cmd_record >/dev/null 2>&1
  GOLDEN="$a" cmd_cli record >/dev/null 2>&1
  GOLDEN="$b" cmd_record >/dev/null 2>&1
  GOLDEN="$b" cmd_cli record >/dev/null 2>&1
  local unstable
  unstable="$(diff -rq "$a" "$b" 2>&1 || true)"
  if [[ -z "$unstable" ]]; then
    ok "goldens are reproducible across two recordings"
  else
    bad "goldens are NOT reproducible -- these contain run-varying data:"
    sed 's/^/          /' <<<"$unstable" | head -8
    while read -r _ f1 _ _; do
      [[ -f "$f1" ]] || continue
      diff "$f1" "${f1/$a/$b}" 2>/dev/null | grep -E '^[<>]' | head -2 | sed 's/^/            /'
    done <<<"$unstable"
  fi
}
# ---------------------------------------------------------------------------

case "${1:-all}" in
  env)     cmd_env ;;
  seeds)   cmd_seeds ;;
  cli)     cmd_cli "${2:-record}" ;;
  record)  cmd_record ;;
  verify)  cmd_verify ;;
  dropped) cmd_dropped ;;
  stable)  cmd_stable ;;
  all)     cmd_env; cmd_seeds; cmd_cli record; cmd_record; cmd_stable; cmd_verify; cmd_dropped ;;
  *) echo "usage: $0 {env|seeds|cli [record|verify]|record|verify|dropped|stable|all}" >&2; exit 2 ;;
esac

echo
echo "==> $PASS passed, $FAIL failed"
[[ $FAIL -eq 0 ]] || exit 1
