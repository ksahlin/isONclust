#!/usr/bin/env bash
#
# Step 0 of repo slimming: get the data out of git before it is stripped.
#
# ---------------------------------------------------------------------------
# This differs from isONcorrect's version of this script, and the difference
# matters. isONcorrect's archived `data/` from the WORKING TREE. isONclust's
# test data was already deleted from HEAD (commits 004e74f, e2be3aa, 2c08b15,
# 710df34, 241d894 and their siblings), so there is nothing in the working tree
# to archive: the only remaining copy of those 519 MB is inside .git, reachable
# through old commits. Once history is rewritten they are gone for good.
#
# So this script extracts blobs from HISTORY, not from disk. It walks every path
# ever committed under the stripped prefixes, pulls the newest blob for each,
# and writes it out with a checksummed manifest.
# ---------------------------------------------------------------------------
#
# By DEFAULT this only builds archives locally and prints what it would upload.
# Uploading publishes ~500 MB under your account and is not something this
# script does on its own -- pass --upload once you have looked at what was built.
#
#   tools/repo-slim/archive_data.sh                 # build locally, report
#   tools/repo-slim/archive_data.sh --upload        # build, then upload to a Release
#
# Options:
#   --outdir DIR   where to build archives (default: ./repo-slim-archive)
#   --tag TAG      release tag (default: data-archive-YYYYMMDD)
#   --upload       actually create the Release and upload
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUTDIR="$REPO_ROOT/repo-slim-archive"
TAG="data-archive-$(date +%Y%m%d)"
UPLOAD=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2 ;;
    --tag)    TAG="$2"; shift 2 ;;
    --upload) UPLOAD=1; shift ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

cd "$REPO_ROOT"

# Everything ever committed under test/. This is the set the rewrite destroys.
# Deliberately NOT read from removal-paths.txt: that file also lists build
# artifacts (__pycache__, .DS_Store) which are regenerable and not worth
# archiving, and analyze.sh has not necessarily been run yet.
# bash 3.2 (what macOS ships) has no `mapfile`, so read into an array the long way.
DATA_PATHS=()
while IFS= read -r line; do
  DATA_PATHS+=("$line")
done < <(
  git log --all --pretty=format: --name-only --no-renames \
    | sed '/^$/d' | LC_ALL=C sort -u | grep '^test/'
)

if [[ ${#DATA_PATHS[@]} -eq 0 ]]; then
  echo "error: no paths under test/ found anywhere in history." >&2
  echo "       Has this repository already been slimmed?" >&2
  exit 1
fi

echo "==> ${#DATA_PATHS[@]} data paths found in history"

RAW="$OUTDIR/raw"
mkdir -p "$RAW/test"

# For each path, take the blob from the most recent commit that HAD it (not the
# commit that deleted it). `git log --all -- <path>` is newest-first; the first
# commit whose tree resolves the path is the last version that existed.
echo "==> extracting blobs from history into $RAW"
TOTAL=0
EXTRACTED=()
for p in "${DATA_PATHS[@]}"; do
  blob=""
  while read -r c; do
    if git cat-file -e "$c:$p" 2>/dev/null; then
      blob="$(git rev-parse "$c:$p")"
      break
    fi
  done < <(git log --all --format=%H -- "$p")

  if [[ -z "$blob" ]]; then
    echo "    WARNING: could not resolve a blob for $p -- skipping" >&2
    continue
  fi

  mkdir -p "$RAW/$(dirname "$p")"
  git cat-file blob "$blob" > "$RAW/$p"
  sz="$(git cat-file -s "$blob")"
  TOTAL=$((TOTAL + sz))
  EXTRACTED+=("$p")
  printf '    %10s bytes  %s  %s\n' "$sz" "${blob:0:8}" "$p"
done

echo "==> extracted ${#EXTRACTED[@]} files, $(python3 -c "print(f'{$TOTAL/1e6:.1f}')") MB"

# ccs.fastq.gz was committed as four `split` parts. Reassembling it is the whole
# point of archiving it -- four opaque parts are not a usable artifact -- so do
# that here and verify it is a valid gzip before throwing the parts away.
if ls "$RAW"/test/ccs.fastq.gz.part-* >/dev/null 2>&1; then
  echo "==> reassembling ccs.fastq.gz from split parts"
  cat "$RAW"/test/ccs.fastq.gz.part-* > "$RAW/test/ccs.fastq.gz"
  if gzip -t "$RAW/test/ccs.fastq.gz" 2>/dev/null; then
    echo "    gzip -t: OK ($(du -h "$RAW/test/ccs.fastq.gz" | cut -f1))"
    rm -f "$RAW"/test/ccs.fastq.gz.part-*
  else
    echo "    gzip -t: FAILED -- keeping the parts and the concatenation" >&2
    echo "    (the parts may have been committed out of order, or truncated)" >&2
  fi
fi

echo "==> packing"
DATA_TAR="$OUTDIR/isONclust-test-data.tar"
tar -cf "$DATA_TAR" -C "$RAW" test
if command -v pigz >/dev/null 2>&1; then pigz -f "$DATA_TAR"; else gzip -f "$DATA_TAR"; fi
DATA_TGZ="${DATA_TAR}.gz"

echo "==> writing manifest"
MANIFEST="$OUTDIR/MANIFEST.txt"
sha() { if command -v shasum >/dev/null 2>&1; then shasum -a 256 "$1"; else sha256sum "$1"; fi; }
{
  echo "isONclust test-data archive"
  echo "Created: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "From commit: $(git rev-parse HEAD)"
  echo
  echo "These files were committed under test/ and are being stripped from git"
  echo "history to make the repository cloneable. They had ALREADY been deleted"
  echo "from HEAD before this archive was made, so at the time of writing the"
  echo "only copy was inside .git, reachable through old commits."
  echo
  echo "## Contents"
  echo "  $(basename "$DATA_TGZ")  $(du -h "$DATA_TGZ" | cut -f1)"
  echo
  echo "## SHA256 of the archive"
  sha "$DATA_TGZ" | sed "s#$OUTDIR/##"
  echo
  echo "## SHA256 of each extracted file"
  ( cd "$RAW" && find test -type f | LC_ALL=C sort | while read -r f; do sha "$f"; done )
  echo
  echo "## Original paths in history, with the blob archived"
  for p in ${EXTRACTED[@]+"${EXTRACTED[@]}"}; do echo "  $p"; done
} > "$MANIFEST"

echo
cat "$MANIFEST"
echo

ASSETS=("$DATA_TGZ" "$MANIFEST")

if [[ $UPLOAD -eq 0 ]]; then
  cat <<EOF
================================================================================
Archives built locally. NOTHING HAS BEEN UPLOADED.

  raw files : $RAW  ($(du -sh "$RAW" | cut -f1))
  archive   : $DATA_TGZ  ($(du -h "$DATA_TGZ" | cut -f1))

Review them, then pick a home.

RECOMMENDED: Zenodo, if any of this data is worth citing.

  1. https://zenodo.org  ->  New upload
  2. Drag in $(basename "$DATA_TGZ") and MANIFEST.txt
  3. Upload type "Dataset"; title e.g.
       "isONclust: test data removed from git history"
  4. Link it to the publication DOI (10.1089/cmb.2019.0299) under
     "Related/alternate identifiers"
  5. Publish, then note the record id from the URL

ALTERNATIVE: a GitHub Release, beside the code, no DOI. Re-run with --upload, or:

  gh release create "$TAG" \\
    --title "Test data removed from git history" \\
    --notes "Data stripped from history to make the repository cloneable. See MANIFEST.txt." \\
$(printf '    %s \\\n' "${ASSETS[@]}" | sed '$ s/ \\$//')

A THIRD OPTION, and the one that was taken: do not publish it at all. None of
these files is referenced by any test, script or README, and the repository's
fixture was replaced by a 356 KB corpus with better discriminating power
(PORTING.md, Finding 5). A local copy is enough.

Note what publishing would NOT achieve: the data is still present in all 9
forks, whose default branches carry it live, and GitHub shares object storage
across a fork network. Stripping it from the canonical repository does not make
it unreachable. See PORTING.md, "Repo hygiene".

Either way: confirm what you keep is readable BEFORE running slim.sh. After the
rewrite, $RAW is the only remaining copy.
================================================================================
EOF
  exit 0
fi

if ! command -v gh >/dev/null 2>&1; then
  echo "error: gh CLI not found; cannot upload." >&2
  exit 1
fi

echo "==> creating release $TAG and uploading $(du -sh "$OUTDIR" | cut -f1)"
gh release create "$TAG" \
  --title "Test data removed from git history" \
  --notes "Data stripped from git history to make the repository cloneable. See MANIFEST.txt for contents and checksums." \
  "${ASSETS[@]}"

echo
echo "==> done. Verify the assets are downloadable BEFORE running slim.sh:"
echo "    gh release view $TAG"
