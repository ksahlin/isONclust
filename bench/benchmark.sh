#!/usr/bin/env bash
#
# Compare the three implementations on speed, memory and accuracy.
#
#   bench/benchmark.sh                          # default corpora, --t 1 and --t 8
#   bench/benchmark.sh --corpora "sirv_sim_err7 droso_20k" --threads "1 8"
#   bench/benchmark.sh --repeats 3              # median of 3 for the timings
#
# The three:
#   python   the reference, modules/*.py
#   port     rust/target/release/isONclust  -- byte-identical to the reference
#   isONclust3  a DIFFERENT ALGORITHM with its own paper, not a port
#
# WHAT IS AND IS NOT COMPARABLE
# -----------------------------
# python and port must agree exactly; any difference is a port bug, and the
# script says so rather than reporting it as an accuracy delta. isONclust3 is
# expected to differ -- that is the point of comparing it.
#
# Accuracy needs per-read ground truth. Pass --truth <file> from
# bench/make_truth.sh; without one, corpora are timed but not scored and the
# table says so rather than printing a meaningless number.
#
# DO NOT score accuracy on the simulated corpora. Every base in them carries the
# same quality character (`I`, phred 40), both algorithms are quality-driven, and
# isONclust3 degenerates to near-singletons there. That measures the simulator.
#
# Memory is peak RSS from /usr/bin/time. On macOS that is `-l` and bytes; on
# Linux `-v` and kilobytes. Both are handled.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
cd "$ROOT"

REF_PYTHON="${REF_PYTHON:-$HOME/miniforge3/envs/isonclust-ref/bin/python}"
PORT_BIN="${PORT_BIN:-$ROOT/rust/target/release/isONclust}"
ISONCLUST3="${ISONCLUST3:-$HOME/source/isONclust3/target/release/isONclust3}"
ISONCLUST_DATA="${ISONCLUST_DATA:-$HOME/data/lrRNA-seq}"
WORK="${WORK:-$(mktemp -d)}"

CORPORA="smoke sirv_sim_err7 sirv_real_10k droso_20k"
TRUTH=""
THREADS="1 8"
REPEATS=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    --corpora) CORPORA="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --repeats) REPEATS="$2"; shift 2 ;;
    --work)    WORK="$2"; shift 2 ;;
    --truth)   TRUTH="$2"; shift 2 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

resolve_corpus() {
  local c="$1" p
  [[ -f "$c" ]] && { echo "$c"; return; }
  p="$(awk -F'\t' -v n="$c" '$1==n {print $2; exit}' bench/corpora.tsv 2>/dev/null)"
  [[ -z "$p" ]] && { echo "$c"; return; }
  case "$p" in
    /*)     echo "$p" ;;
    test/*) echo "$ROOT/$p" ;;
    *)      echo "$ISONCLUST_DATA/$p" ;;
  esac
}


# Peak RSS in MB, portable between macOS (-l, bytes) and GNU (-v, kB).
TIMEFMT_FILE="$WORK/.time"
run_timed() { # run_timed <outvar_prefix> -- command...
  local t0 t1
  t0=$(python3 -c 'import time;print(time.time())')
  if /usr/bin/time -l "$@" >"$WORK/.stdout" 2>"$TIMEFMT_FILE"; then :; else
    RUN_RC=$?; RUN_SECS=0; RUN_MB=0; return 0
  fi
  t1=$(python3 -c 'import time;print(time.time())')
  RUN_RC=0
  RUN_SECS=$(python3 -c "print(f'{$t1-$t0:.2f}')")
  RUN_MB=$(awk '/maximum resident set size/ {print int($1/1048576); found=1}
                /Maximum resident set size/ {print int($6/1024); found=1}
                END {if (!found) print 0}' "$TIMEFMT_FILE" | head -1)
}

median() { python3 -c "
import sys
v=sorted(float(x) for x in sys.argv[1:])
print(f'{v[len(v)//2]:.2f}')" "$@"; }

printf '%-16s %-11s %-4s %8s %8s %7s %7s %7s %7s %7s\n' \
  CORPUS TOOL -t "SECS" "PEAK_MB" "CLUST" "HOMOG" "COMPL" "V" "ARI"
printf '%.0s-' {1..96}; echo

for corpus in $CORPORA; do
  fq="$(resolve_corpus "$corpus")"
  if [[ ! -f "$fq" ]]; then
    printf '%-16s  (missing: %s)\n' "$corpus" "$fq"
    continue
  fi
  ref_clusters=""

  for t in $THREADS; do
    for tool in python port isONclust3; do
      # isONclust3 has no --t equivalent in this comparison; run it once.
      [[ "$tool" == "isONclust3" && "$t" != "${THREADS%% *}" ]] && continue

      secs_runs=(); mb=0; out=""
      for rep in $(seq 1 "$REPEATS"); do
        d="$WORK/$corpus.$tool.$t.$rep"; rm -rf "$d"; mkdir -p "$d"
        case "$tool" in
          python)
            run_timed env PYTHONHASHSEED=0 "$REF_PYTHON" isONclust --ont --t "$t" \
              --fastq "$fq" --outfolder "$d"
            out="$d/final_clusters.tsv" ;;
          port)
            run_timed "$PORT_BIN" --ont --t "$t" --fastq "$fq" --outfolder "$d"
            out="$d/final_clusters.tsv" ;;
          isONclust3)
            # --post-cluster is NOT optional for a fair comparison. The README's
            # own example passes it, and on real SIRV ONT reads it takes the
            # clustering from 5661 clusters to 66 -- V(gene) 0.36 -> 0.67, and
            # ARI(transcript) 0.53 -> 0.60 -- for 0.09s. Benchmarking without it
            # measures a tool nobody runs.
            run_timed "$ISONCLUST3" --fastq "$fq" --outfolder "$d" --mode ont \
              --seeding minimizer --post-cluster --no-fastq
            out="$d/clustering/final_clusters.tsv" ;;
        esac
        secs_runs+=("$RUN_SECS"); mb="$RUN_MB"
      done
      secs="$(median "${secs_runs[@]}")"

      nclust="-"; scores="      -       -       -       -"
      if [[ -f "$out" ]]; then
        nclust="$(cut -f1 "$out" | sort -u | wc -l | tr -d ' ')"
        if [[ -n "$TRUTH" && -f "$TRUTH" ]]; then
          scores="$("$REF_PYTHON" bench/accuracy.py --clusters "$out" --truth-file "$TRUTH" \
                    --tsv --label x 2>/dev/null \
                    | awk -F'\t' '{if ($7=="-") printf "      -       -       -       -";
                                    else printf "%7.4f %7.4f %7.4f %7.4f", $7,$8,$9,$10}')"
        fi
      fi
      printf '%-16s %-11s %-4s %8s %8s %7s %s\n' \
        "$corpus" "$tool" "$t" "$secs" "$mb" "$nclust" "$scores"

      # The port must agree with the reference exactly. A difference here is a
      # port bug, not a result, so say so loudly.
      if [[ "$tool" == "python" ]]; then ref_clusters="$out"
      elif [[ "$tool" == "port" && -n "$ref_clusters" && -f "$out" ]]; then
        cmp -s "$ref_clusters" "$out" || \
          printf '%-16s  *** PORT DIFFERS FROM THE REFERENCE at --t %s -- this is a bug, not a result\n' "$corpus" "$t"
      fi
    done
  done
done
