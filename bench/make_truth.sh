#!/usr/bin/env bash
#
# Derive per-read ground truth for a real corpus by aligning it to a reference.
#
#   bench/make_truth.sh READS.fastq REFERENCE.fasta OUT.tsv [--level gene|transcript]
#
# Writes `<read accession>\t<class>` for every read with a primary alignment.
#
# WHY THIS EXISTS
# ---------------
# The simulated SIRV corpora carry their source transcript in the accession, so
# they need no alignment -- but they are USELESS for comparing tools, because
# every base carries the same quality character (`I`, phred 40). Both isONclust
# and isONclust3 are quality-driven algorithms, and on constant quality
# isONclust3 degenerates to near-singletons. Scoring accuracy there measures the
# simulator, not the tool. See PORTING.md.
#
# So cross-tool accuracy has to be scored on REAL reads, and real reads need
# their truth derived by alignment. This is what the paper's own
# `scripts/compute_cluster_quality.py` does, via a sorted BAM; this is the same
# idea with minimap2 and no pysam dependency.
#
# LEVEL
#   gene       SIRV101 -> SIRV1     THE TARGET. isONclust is a gene clustering
#                                   tool -- "each cluster represents all reads
#                                   that came from a gene" -- so this is the
#                                   level any result should be judged at, and it
#                                   is the default.
#   transcript SIRV101 -> SIRV101   DIAGNOSTIC ONLY. Useful for understanding the
#                                   shape of a clustering; never a thing to
#                                   optimise. It disagrees with the gene-level
#                                   verdict on ONT and agrees on PacBio, so
#                                   steering by it would mean steering by a
#                                   number that changes its mind per platform.
#                                   See PORTING.md, "Accuracy".
set -euo pipefail

READS="${1:?reads fastq}"
REF="${2:?reference fasta}"
OUT="${3:?output tsv}"
LEVEL="gene"
shift 3 || true
while [[ $# -gt 0 ]]; do
  case "$1" in
    --level) LEVEL="$2"; shift 2 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

command -v minimap2 >/dev/null || { echo "error: minimap2 not found" >&2; exit 1; }

# -x map-ont for nanopore; --secondary=no so each read has at most one hit.
minimap2 -ax map-ont --secondary=no -t 8 "$REF" "$READS" 2>/dev/null \
  | awk -v level="$LEVEL" '
      /^@/ { next }
      {
        flag = $2
        # Skip unmapped (4), secondary (256) and supplementary (2048).
        # Arithmetic rather than and(): macOS ships BWK awk, which has no
        # bitwise functions, and this script has to run on the machine it is on.
        if (int(flag/4) % 2 || int(flag/256) % 2 || int(flag/2048) % 2) next
        acc = $1; ref = $3
        if (ref == "*") next
        cls = ref
        if (level == "gene") {
          # SIRV101 -> SIRV1; anything else is left alone
          if (match(ref, /^SIRV[0-9]/)) cls = substr(ref, 1, 5)
        }
        print acc "\t" cls
      }' | LC_ALL=C sort -u > "$OUT.qname"

# The tools disagree on what a read is called, so emit truth under every form
# they use. minimap2's QNAME stops at the first whitespace; isONclust's readfq
# replaces spaces with underscores and keeps the whole header; isONclust3 keeps
# the first token. Keying on only one of those silently scores zero reads.
awk 'NR%4==1' "$READS" | sed 's/^@//' > "$OUT.headers"
awk -F'\t' '
  NR==FNR { cls[$1] = $2; next }
  {
    full = $0
    split(full, parts, " ")
    qname = parts[1]
    if (qname in cls) {
      under = full
      gsub(/ /, "_", under)
      print qname "\t" cls[qname]
      if (under != qname) print under "\t" cls[qname]
    }
  }' "$OUT.qname" "$OUT.headers" | LC_ALL=C sort -u > "$OUT"
rm -f "$OUT.qname" "$OUT.headers"

echo "  wrote $OUT: $(wc -l < "$OUT" | tr -d ' ') accession forms covering $(cut -f2 "$OUT" | sort -u | wc -l | tr -d ' ') classes (level: $LEVEL)"
