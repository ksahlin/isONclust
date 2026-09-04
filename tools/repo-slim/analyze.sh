#!/usr/bin/env bash
#
# Step 1 of repo slimming: work out what would be removed, and prove it is worth
# doing, without changing anything.
#
# Writes two files next to this script:
#   removal-paths.txt   the exact path list slim.sh will strip (review this)
#   analysis.txt        a human-readable report
#
# Nothing here mutates the repository. Run it, read removal-paths.txt, and only
# then run slim.sh.
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REMOVAL="$HERE/removal-paths.txt"
REPORT="$HERE/analysis.txt"

# ---------------------------------------------------------------------------
# What survives, and what does not.
#
# isONclust's situation differs from isONcorrect's in one way that changes this
# list: every file under test/ had ALREADY been deleted from HEAD before the
# rewrite (commits 004e74f, e2be3aa, 2c08b15, 710df34, 241d894 and siblings).
# So there is no fixture in the working tree to protect, and KEEP is empty.
#
# That is a deliberate choice, not an oversight. archive_data.sh extracts every
# stripped file to repo-slim-archive/raw/test/ first, so nothing is lost, and
# keeping a historical blob would not have put the file back in HEAD anyway --
# which is where a fixture needs to be.
#
# In the event the old fixture was NOT restored: sample_alz_2k.fastq turned out
# to be blind to the minimizer window (PORTING.md, Finding 5), and was replaced
# by test/sirv_sim_120.fastq (356 KB, simulated SIRV with ground truth) plus real
# corpora kept outside the repository. See bench/corpora.tsv.
KEEP=()

# Prefixes whose contents are stripped from all history unless listed in KEEP.
#
#   test/                  -- 519 MB, the entire reason for this exercise
#   modules/__pycache__/   -- committed CPython bytecode, 1.15 MB across
#                             versions, regenerable, and wrong for any
#                             interpreter but 3.6
STRIP_PREFIXES=("test/" "modules/__pycache__/")

# Paths stripped by basename rather than prefix. macOS Finder metadata was
# committed at two paths (modules/.DS_Store and scripts/.DS_Store, which share
# one blob) and neither is under a strippable prefix.
STRIP_BASENAMES=(".DS_Store")

cd "$REPO_ROOT"

echo "==> enumerating every object in history (a few seconds on a 500 MB repo)"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

git rev-list --objects --all > "$TMP/allobj.txt"
awk 'NF>1' "$TMP/allobj.txt" > "$TMP/withpath.txt"
cut -d' ' -f1 "$TMP/withpath.txt" \
  | git cat-file --batch-check='%(objectname) %(objecttype) %(objectsize)' > "$TMP/sizes.txt"

# The removal list MUST come from a tree walk, not from the object listing.
# `git rev-list --objects` emits each blob exactly once with a single path, so
# any file whose content is byte-identical to another file is invisible there --
# test_data/chr6_ensemble.fa is deduplicated against data/chr6_transcripts.fa and
# would silently survive the rewrite. This enumerates every path ever committed.
echo "==> enumerating every path ever committed"
git log --all --pretty=format: --name-only --no-renames \
  | sed '/^$/d' | LC_ALL=C sort -u > "$TMP/allpaths.txt"
echo "    $(wc -l < "$TMP/allpaths.txt" | tr -d ' ') distinct paths in history"

echo "==> classifying"
KEEP_JOINED="$(printf '%s\n' ${KEEP[@]+"${KEEP[@]}"})"
PREFIX_JOINED="$(printf '%s\n' "${STRIP_PREFIXES[@]}")"

BASENAME_JOINED="$(printf '%s\n' "${STRIP_BASENAMES[@]}")"

KEEP_LIST="$KEEP_JOINED" PREFIX_LIST="$PREFIX_JOINED" BASENAME_LIST="$BASENAME_JOINED" \
python3 - "$TMP/withpath.txt" "$TMP/sizes.txt" "$REMOVAL" "$REPORT" "$TMP/allpaths.txt" <<'PY'
import os, sys
from collections import defaultdict

withpath, sizesf, removal_out, report_out, allpathsf = sys.argv[1:6]
keep = {l for l in os.environ["KEEP_LIST"].splitlines() if l}
prefixes = tuple(l for l in os.environ["PREFIX_LIST"].splitlines() if l)
basenames = tuple(l for l in os.environ["BASENAME_LIST"].splitlines() if l)

sizes = {}
for line in open(sizesf):
    p = line.split()
    if len(p) == 3 and p[1] == "blob":
        sizes[p[0]] = int(p[2])

# A blob can be reachable at several paths and git stores it once, so attribute
# each distinct blob to the first path that reaches it and never double count.
# The authoritative removal list: every path ever committed under a stripped
# prefix that is not explicitly kept. Derived from the tree walk, so
# content-deduplicated paths are included.
strip_paths = {
    p for p in (l.rstrip("\n") for l in open(allpathsf))
    if p and (p.startswith(prefixes) or p.rsplit("/", 1)[-1] in basenames)
    and p not in keep
}

seen = set()
strip_bytes = defaultdict(int)
keep_bytes = defaultdict(int)
other_total = 0
total = 0

# `git rev-list --objects` emits each object exactly once, with whichever path
# reached it first, so this loop cannot be used to decide whether a given path
# exists -- deduplicated paths simply never appear. Fixture existence is checked
# separately, against the trees, after this step.
for line in open(withpath):
    sha, _, path = line.partition(" ")
    path = path.strip()
    if sha not in sizes or sha in seen:
        continue
    seen.add(sha)
    size = sizes[sha]
    total += size
    if path in strip_paths:
        strip_bytes[path] += size
    elif path in keep:
        keep_bytes[path] += size
    else:
        other_total += size

with open(removal_out, "w") as fh:
    fh.write("# Paths stripped from ALL history by tools/repo-slim/slim.sh.\n")
    fh.write("# Generated by analyze.sh -- review before running slim.sh.\n")
    fh.write("# Format is git-filter-repo --paths-from-file: one literal path per line.\n")
    for p in sorted(strip_paths):
        fh.write(p + "\n")

stripped_total = sum(strip_bytes.values())
kept_total = sum(keep_bytes.values())

lines = []
w = lines.append
w("Repository slimming analysis")
w("=" * 60)
w("")
w(f"total distinct blob bytes in history : {total/1e9:8.3f} GB")
w(f"  to be stripped                     : {stripped_total/1e9:8.3f} GB  ({100*stripped_total/total:.1f}%)")
w(f"  retained (source, scripts, docs)   : {other_total/1e9:8.3f} GB")
w(f"  retained fixtures (KEEP list)      : {kept_total/1e6:8.3f} MB")
w("")
w(f"paths stripped: {len(strip_paths)}")
w("")
w("Largest paths being removed")
w("-" * 60)
for p, s in sorted(strip_bytes.items(), key=lambda x: -x[1])[:25]:
    w(f"  {s/1e6:10.2f} MB  {p}")
w("")
w("Note: git stores each distinct blob once, so identical files at different")
w("paths are counted under whichever path reaches them first. This is not")
w("cosmetic -- it reflects real storage. Checked on this history:")
w("")
w("  * Only three blobs in all of isONclust's history are reachable at more")
w("    than one path: modules/.DS_Store == scripts/.DS_Store (both stripped),")
w("    cemetary/cluster_parallel.py == modules/cluster_parallel.py, and")
w("    modules/compute_shared_minimizers_probabilities.py ==")
w("    scripts/compute_shared_minimizer_probabilities.py (all retained). So the")
w("    deduplication trap that made isONcorrect's first removal list wrong does")
w("    not bite here -- but this list is still built from a tree walk, because")
w("    that is the only way to know that.")
w("")
w("  * test/ccs.fastq.gz was committed as four `split` parts, part-aa..part-ad,")
w("    365 MB together. They are four blobs, not one file, and archive_data.sh")
w("    reassembles them (verified with `gzip -t`) before they are destroyed.")
w("")
w("  * KEEP is empty on purpose. Everything under test/ was already deleted")
w("    from HEAD, so protecting a historical blob would not restore it. The")
w("    replacement fixture is committed fresh after the rewrite; see")
w("    bench/corpora.tsv.")

report = "\n".join(lines)
open(report_out, "w").write(report + "\n")
print(report)
PY

# ---------------------------------------------------------------------------
# Fixture verification, done against the trees rather than the object listing.
#
# A KEEP entry naming a path that does not exist protects nothing and would let
# a fixture be stripped silently, so this is a hard error.
# ---------------------------------------------------------------------------
{
  echo
  echo "Fixtures explicitly retained"
  echo "------------------------------------------------------------"
  [[ ${#KEEP[@]} -eq 0 ]] && echo "  (none -- see the KEEP comment in analyze.sh)"
} | tee -a "$REPORT"

FIXTURE_ERROR=0
for p in ${KEEP[@]+"${KEEP[@]}"}; do
  if git cat-file -e "HEAD:$p" 2>/dev/null; then
    blob="$(git rev-parse "HEAD:$p")"
    size="$(git cat-file -s "$blob")"
    # Other paths at HEAD sharing this exact blob. `|| true` because grep exits
    # non-zero when a fixture has no twins, which is the normal case.
    twins="$(git ls-tree -r HEAD | awk -v b="$blob" '$3==b {print $4}' | { grep -vFx "$p" || true; } | paste -sd', ' -)"
    line="$(printf '  %10.2f KB  %s' "$(echo "$size/1000" | bc -l)" "$p")"
    [[ -n "$twins" ]] && line="$line"$'\n'"             shares content with: $twins"
  elif [[ -n "$(git log --all --format=%H -1 -- "$p" 2>/dev/null)" ]]; then
    line="  $(printf '%10s' 'history')  $p   (not in HEAD, but present in history)"
  else
    line="  $(printf '%10s' 'MISSING')  $p   <-- ERROR: path not found in HEAD or history"
    FIXTURE_ERROR=1
  fi
  echo "$line" | tee -a "$REPORT"
done

if [[ $FIXTURE_ERROR -eq 1 ]]; then
  {
    echo
    echo "ERROR: the KEEP list names paths that do not exist. Fix the KEEP array in"
    echo "analyze.sh before running slim.sh -- those entries protect nothing."
  } | tee -a "$REPORT"
fi

echo
echo "==> wrote $REMOVAL"
echo "==> wrote $REPORT"
echo

if [[ $FIXTURE_ERROR -eq 1 ]]; then
  echo "REFUSING to proceed: fix the KEEP list first." >&2
  exit 1
fi

echo "Review removal-paths.txt, then run: tools/repo-slim/slim.sh"
