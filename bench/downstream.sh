#!/usr/bin/env bash
#
# The ultimate quality check: does a better clustering give better transcripts?
#
#   bench/downstream.sh --corpus sirv_pacbio --out DIR
#
# Every metric before this one is intrinsic -- it scores the clustering against a
# truth set. This scores what the clustering is FOR: the isoforms isONform
# reconstructs downstream. A clustering can look good by V-measure and still hand
# isONform the wrong reads.
#
# PIPELINE. On PacBio the reads are accurate enough to go straight from
# clustering to isoform reconstruction:
#
#     isONclust  ->  isONform
#
# On ONT they are not, and isONcorrect has to run in between:
#
#     isONclust  ->  isONcorrect  ->  isONform
#
# which is why this script is run on PacBio: one fewer heavy stage, and it is the
# platform where the intrinsic metrics already separate the two clusterers most
# clearly.
#
# WHAT IS COMPARED. isONclust1 and this port produce byte-identical clusterings,
# so the downstream is run TWICE, not three times: isONclust1 (standing for both)
# against isONclust3. Running the port separately would measure nothing.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
cd "$ROOT"

PORT_BIN="${PORT_BIN:-$ROOT/rust/target/release/isONclust}"
ISONCLUST3="${ISONCLUST3:-$HOME/source/isONclust3/target/release/isONclust3}"
ISONFORM="${ISONFORM:-$HOME/source/isONform/rust/target/release/isONform_parallel}"
ISONFORM_SCORER="${ISONFORM_SCORER:-$HOME/source/isONform/bench/accuracy_isoforms.py}"
REF_PYTHON="${REF_PYTHON:-$HOME/miniforge3/envs/isonclust-ref/bin/python}"
ISONCLUST_DATA="${ISONCLUST_DATA:-$HOME/data/lrRNA-seq}"

CORPUS="sirv_pacbio"
OUT=""
MINREADS=3
THREADS=8
while [[ $# -gt 0 ]]; do
  case "$1" in
    --corpus)   CORPUS="$2"; shift 2 ;;
    --out)      OUT="$2"; shift 2 ;;
    --min-reads) MINREADS="$2"; shift 2 ;;
    --t)        THREADS="$2"; shift 2 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done
[[ -n "$OUT" ]] || OUT="$(mktemp -d)"
mkdir -p "$OUT"

resolve() {
  local p; p="$(awk -F'\t' -v n="$1" '$1==n {print $2; exit}' bench/corpora.tsv)"
  case "$p" in /*) echo "$p" ;; test/*) echo "$ROOT/$p" ;; *) echo "$ISONCLUST_DATA/$p" ;; esac
}
preset() { awk -F'\t' -v n="$1" '$1==n {print $4; exit}' bench/corpora.tsv; }

FQ="$(resolve "$CORPUS")"
PRESET="$(preset "$CORPUS")"
case "$PRESET" in isoseq) I1="--isoseq"; I3MODE="pacbio" ;; *) I1="--ont"; I3MODE="ont" ;; esac
echo "==> $CORPUS ($(( $(wc -l < "$FQ") / 4 )) reads), preset $PRESET"

# --- 1. cluster -------------------------------------------------------------
echo "==> clustering"
c1="$OUT/isonclust1"; rm -rf "$c1"; mkdir -p "$c1"
"$PORT_BIN" $I1 --t "$THREADS" --fastq "$FQ" --outfolder "$c1" >/dev/null 2>&1
echo "    isONclust1: $(cut -f1 "$c1/final_clusters.tsv" | sort -u | wc -l | tr -d ' ') clusters"

c3="$OUT/isonclust3"; rm -rf "$c3"; mkdir -p "$c3"
"$ISONCLUST3" --fastq "$FQ" --outfolder "$c3" --mode "$I3MODE" \
  --seeding minimizer --post-cluster -n "$MINREADS" >/dev/null 2>&1
echo "    isONclust3: $(cut -f1 "$c3/clustering/final_clusters.tsv" | sort -u | wc -l | tr -d ' ') clusters"

# --- 2. per-cluster fastq ---------------------------------------------------
# isONform consumes a folder of per-cluster fastq files. isONclust3 writes them
# itself; isONclust1 needs its write_fastq subcommand. Both are filtered to
# clusters of at least --min-reads, because a 1-read cluster cannot produce a
# consensus and would only add noise to both sides equally.
echo "==> writing per-cluster fastq (clusters with >= $MINREADS reads)"
f1="$OUT/fastq1"; rm -rf "$f1"
"$PORT_BIN" write_fastq --clusters "$c1/final_clusters.tsv" --fastq "$FQ" \
  --outfolder "$f1" --N "$MINREADS" >/dev/null 2>&1
echo "    isONclust1: $(ls "$f1" | wc -l | tr -d ' ') cluster files"

f3="$OUT/fastq3"; rm -rf "$f3"; mkdir -p "$f3"
if [[ -d "$c3/clustering/fastq_files" ]]; then
  for f in "$c3/clustering/fastq_files"/*.fastq; do
    [[ -f "$f" ]] || continue
    if [[ $(( $(wc -l < "$f") / 4 )) -ge $MINREADS ]]; then cp "$f" "$f3/"; fi
  done
fi
echo "    isONclust3: $(ls "$f3" | wc -l | tr -d ' ') cluster files"

# --- 3. isoform reconstruction ----------------------------------------------
for tag in 1 3; do
  src="$OUT/fastq$tag"; dst="$OUT/isoforms$tag"; rm -rf "$dst"; mkdir -p "$dst"
  n=$(ls "$src" 2>/dev/null | wc -l | tr -d ' ')
  [[ "$n" -gt 0 ]] || { echo "==> isONform on isONclust$tag: no clusters, skipping"; continue; }
  echo "==> isONform on isONclust$tag ($n clusters)"
  t0=$(python3 -c 'import time;print(time.time())')
  "$ISONFORM" --fastq_folder "$src" --outfolder "$dst" --t "$THREADS" >"$dst.log" 2>&1 || true
  t1=$(python3 -c 'import time;print(time.time())')
  python3 -c "print(f'    {$t1-$t0:.1f}s')"
  find "$dst" -name '*.fa*' | head -3 | sed 's/^/      /'
done

# --- 4. score the reconstructed isoforms ------------------------------------
#
# isONform writes one `cluster<N>_merged.fa` per cluster; the reconstructed
# transcriptome is their concatenation. Names are made unique by prefixing the
# cluster, because the same isoform id recurs across clusters.
echo "==> collecting isoforms"
for tag in 1 3; do
  dst="$OUT/isoforms$tag"; out="$OUT/transcriptome$tag.fa"; : > "$out"
  n=0
  for f in "$dst"/cluster*_merged.fa; do
    [[ -f "$f" ]] || continue
    cl="$(basename "$f" | sed 's/_merged\.fa$//')"
    awk -v c="$cl" '/^>/ {print ">" c "_" substr($0,2); next} {print}' "$f" >> "$out"
    n=$((n+1))
  done
  echo "    isONclust$tag: $(grep -c '^>' "$out" 2>/dev/null || echo 0) isoforms from $n clusters -> $out"
done

REF_TRANSCRIPTOME="${REF_TRANSCRIPTOME:-$ISONCLUST_DATA/sirv/sirv_transcriptome.fasta}"
if [[ -f "$ISONFORM_SCORER" && -f "$REF_TRANSCRIPTOME" ]]; then
  echo "==> scoring against $(basename "$REF_TRANSCRIPTOME")"
  "$REF_PYTHON" "$ISONFORM_SCORER" \
    --transcriptome "$REF_TRANSCRIPTOME" \
    --isoforms "isONclust1=$OUT/transcriptome1.fa" \
    --isoforms "isONclust3=$OUT/transcriptome3.fa" 2>&1 | sed 's/^/    /'
else
  echo "==> scorer or reference transcriptome missing; isoforms are in $OUT"
fi

echo
echo "==> outputs are under $OUT"
