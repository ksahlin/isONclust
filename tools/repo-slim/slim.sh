#!/usr/bin/env bash
#
# Step 2 of repo slimming: rewrite history in a throwaway clone and verify it.
#
# This script NEVER touches your working repository and NEVER pushes. It produces
# a rewritten mirror in a scratch directory, checks it, and prints the commands
# you would run to publish it. Pushing is your call and yours alone -- it is the
# irreversible part.
#
#   tools/repo-slim/analyze.sh      # first: produce and review removal-paths.txt
#   tools/repo-slim/slim.sh         # then: rewrite + verify
#
# Options:
#   --workdir DIR   where to build the rewritten mirror (default: a temp dir)
#   --keep-workdir  don't delete the workdir on failure
#
# Requires git-filter-repo:
#   pipx install git-filter-repo      (or)   pip install git-filter-repo
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REMOVAL="$HERE/removal-paths.txt"

WORKDIR=""
KEEP_WORKDIR=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --workdir)      WORKDIR="$2"; shift 2 ;;
    --keep-workdir) KEEP_WORKDIR=1; shift ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

# --- preconditions -----------------------------------------------------------

if [[ ! -f "$REMOVAL" ]]; then
  echo "error: $REMOVAL not found. Run tools/repo-slim/analyze.sh first." >&2
  exit 1
fi

FILTER_REPO="${GIT_FILTER_REPO:-}"
if [[ -z "$FILTER_REPO" ]]; then
  if command -v git-filter-repo >/dev/null 2>&1; then
    FILTER_REPO="$(command -v git-filter-repo)"
  else
    echo "error: git-filter-repo not found." >&2
    echo "  install with:  pipx install git-filter-repo" >&2
    echo "            or:  pip install git-filter-repo" >&2
    echo "  or set GIT_FILTER_REPO=/path/to/git-filter-repo" >&2
    exit 1
  fi
fi

ORIGIN="$(git -C "$REPO_ROOT" remote get-url origin)"
echo "==> source repo:  $REPO_ROOT"
echo "==> origin:       $ORIGIN"
echo "==> filter-repo:  $FILTER_REPO"
echo "==> removing:     $(grep -vc '^#' "$REMOVAL") paths"
echo

# --- build a throwaway mirror ------------------------------------------------

CLEANUP_WORKDIR=0
if [[ -z "$WORKDIR" ]]; then
  WORKDIR="$(mktemp -d)"
  CLEANUP_WORKDIR=1
fi
mkdir -p "$WORKDIR"
MIRROR="$WORKDIR/isONclust-slim.git"

cleanup() {
  if [[ $? -ne 0 && $KEEP_WORKDIR -eq 0 && $CLEANUP_WORKDIR -eq 1 ]]; then
    rm -rf "$WORKDIR"
  fi
}
trap cleanup EXIT

echo "==> cloning a fresh mirror into $MIRROR"
echo "    (filter-repo requires a pristine clone; your working repo is untouched)"
rm -rf "$MIRROR"
git clone --mirror --no-local "$REPO_ROOT" "$MIRROR" 2>&1 | sed 's/^/    /'

size_of() { du -sk "$1" | awk '{printf "%.0f", $1/1024}'; }
BEFORE_MB="$(size_of "$MIRROR")"
echo "==> size before: ${BEFORE_MB} MB"

# --- rewrite -----------------------------------------------------------------

echo
echo "==> rewriting history"
# --invert-paths: the file lists what to REMOVE.
# --force: the mirror is freshly made by us, so filter-repo's freshness guard is
#          satisfied but the mirror carries origin refs it wants confirmation on.
(
  cd "$MIRROR"
  "$FILTER_REPO" --paths-from-file "$REMOVAL" --invert-paths --force 2>&1 | sed 's/^/    /'
  git reflog expire --expire=now --all
  git gc --prune=now --aggressive --quiet
)

AFTER_MB="$(size_of "$MIRROR")"
echo "==> size after:  ${AFTER_MB} MB"

# --- verify ------------------------------------------------------------------

echo
echo "==> verifying"
FAIL=0
check() { # check <description> <expected> <actual>
  if [[ "$2" == "$3" ]]; then
    printf '    ok    %-52s %s\n' "$1" "$3"
  else
    printf '    FAIL  %-52s expected %s, got %s\n' "$1" "$2" "$3"
    FAIL=1
  fi
}

# blob_at <repo> <path> -> blob sha, or the literal "absent". `git rev-parse`
# prints the unresolved string on stdout when a path is missing, so it cannot be
# used here; cat-file -e is the reliable existence test.
blob_at() {
  if git -C "$1" cat-file -e "HEAD:$2" 2>/dev/null; then
    git -C "$1" rev-parse "HEAD:$2"
  else
    echo absent
  fi
}

# filter-repo drops commits that become empty once their only content is
# stripped, so the total count legitimately falls. What must NOT change is the
# number of commits touching source.
ORIG_COMMITS="$(git -C "$REPO_ROOT" rev-list --count --all)"
NEW_COMMITS="$(git -C "$MIRROR" rev-list --count --all)"
printf '    info  %-52s %s -> %s (%s data-only commits pruned)\n' \
  "commit count" "$ORIG_COMMITS" "$NEW_COMMITS" "$((ORIG_COMMITS - NEW_COMMITS))"

SRC_PATHSPEC=(isONclust modules scripts cemetary README.md setup.py setup.cfg requirements.txt MANIFEST.in LICENSE.txt .travis.yml)
ORIG_SRC="$(git -C "$REPO_ROOT" rev-list --count --all -- "${SRC_PATHSPEC[@]}")"
NEW_SRC="$(git -C "$MIRROR" rev-list --count --all -- "${SRC_PATHSPEC[@]}")"
check "commits touching source preserved" "$ORIG_SRC" "$NEW_SRC"

# Source must survive untouched. Every .py in the tool, plus the packaging and
# docs. p_minimizers_shared.py is included deliberately: it is a 2.5 MB
# generated Python literal and the largest retained file, so it is the one most
# likely to be swept up by a careless prefix rule.
for f in isONclust \
         modules/__init__.py modules/cluster.py modules/consensus.py \
         modules/get_sorted_fastq_for_cluster.py modules/help_functions.py \
         modules/p_minimizers_shared.py modules/parallelize.py \
         scripts/compute_cluster_quality.py \
         scripts/compute_shared_minimizer_probabilities.py \
         cemetary/cluster_parallel.py cemetary/old_functions.py \
         setup.py setup.cfg requirements.txt MANIFEST.in LICENSE.txt \
         README.md .travis.yml; do
  check "unchanged: $f" "$(blob_at "$REPO_ROOT" "$f")" "$(blob_at "$MIRROR" "$f")"
done

# No fixtures are retained -- KEEP is empty (see analyze.sh). The one file that
# is replaced (test/sirv_sim_120.fastq) is committed fresh after the rewrite;
# archive_data.sh has already extracted everything. Assert the old paths really
# are gone rather than assuming it, so a silent partial strip cannot pass.

# Stripped paths must be gone from HEAD. The test/ entries were already deleted
# from HEAD before the rewrite, so those assertions are vacuous today and are
# kept only so this script stays correct if it is ever re-run on a tree that
# still has them. The __pycache__ and .DS_Store entries are NOT vacuous: they
# are in HEAD right now and the rewrite is what removes them.
for f in modules/.DS_Store scripts/.DS_Store \
         modules/__pycache__/cluster.cpython-36.pyc \
         modules/__pycache__/p_minimizers_shared.cpython-36.pyc \
         test/sample_alz_2k.fastq test/ccs.fastq.gz.part-aa; do
  check "stripped from HEAD: $f" "absent" "$(blob_at "$MIRROR" "$f")"
done

# ...and from all of history. Walk the trees, not the object listing, for the
# same dedup reason that made the removal list wrong when built from objects.
# KEEP is empty, so nothing under the stripped prefixes may survive at all and
# there is no keep-list regex to subtract.
LEFTOVER="$(git -C "$MIRROR" log --all --pretty=format: --name-only --no-renames \
            | sed '/^$/d' | LC_ALL=C sort -u \
            | grep -E '^(test/|modules/__pycache__/)|(^|/)\.DS_Store$' || true)"
if [[ -z "$LEFTOVER" ]]; then
  printf '    ok    %-52s %s\n' "no stripped paths anywhere in history" "clean"
else
  printf '    FAIL  %-52s %s remain\n' "stripped paths still in history" "$(wc -l <<<"$LEFTOVER" | tr -d ' ')"
  head -5 <<<"$LEFTOVER" | sed 's/^/          /'
  FAIL=1
fi

# The 0.0.4 tag's tree contains test/sample_alz_2k.fastq, both .DS_Store files
# and all four .pyc files, so filter-repo has to rewrite it. Two ways this can
# go wrong and both are silent: the tag can vanish, or it can survive pointing
# at an unrewritten commit -- which would keep every stripped byte reachable on
# the server AND get pulled by `git clone`, which fetches tags by default,
# undoing the entire exercise.
ORIG_TAGS="$(git -C "$REPO_ROOT" tag -l | LC_ALL=C sort | paste -sd, -)"
NEW_TAGS="$(git -C "$MIRROR" tag -l | LC_ALL=C sort | paste -sd, -)"
check "tags preserved" "$ORIG_TAGS" "$NEW_TAGS"

for t in $(git -C "$MIRROR" tag -l); do
  DIRTY="$(git -C "$MIRROR" ls-tree -r --name-only "$t^{tree}" \
           | grep -E '^(test/|modules/__pycache__/)|(^|/)\.DS_Store$' || true)"
  if [[ -z "$DIRTY" ]]; then
    printf '    ok    %-52s %s\n' "tag tree clean: $t" "clean"
  else
    printf '    FAIL  %-52s %s stripped paths remain\n' "tag tree dirty: $t" "$(wc -l <<<"$DIRTY" | tr -d ' ')"
    sed 's/^/          /' <<<"$DIRTY"
    FAIL=1
  fi
done

# The whole point.
SHRINK="$(python3 -c "print(f'{$BEFORE_MB/max($AFTER_MB,1):.0f}x')")"
printf '    info  %-52s %s MB -> %s MB (%s smaller)\n' "repository size" "$BEFORE_MB" "$AFTER_MB" "$SHRINK"

echo
if [[ $FAIL -ne 0 ]]; then
  echo "VERIFICATION FAILED — do not push this. Workdir kept at:" >&2
  echo "  $MIRROR" >&2
  KEEP_WORKDIR=1
  exit 1
fi
echo "==> verification passed"

# --- hand over ---------------------------------------------------------------

cat <<EOF

================================================================================
Rewritten mirror is ready. NOTHING HAS BEEN PUSHED.

  $MIRROR

Before you push, understand what this does:

  * Every commit SHA changes. All 9 forks diverge permanently and cannot be
    fast-forwarded. Anyone with a clone (80 stargazers, 9 forks) must re-clone.
  * The 0.0.4 tag moves to a rewritten commit. Verified above that it survives
    and that its tree is clean.
  * Commit-pinned links break, including any in the two papers
    (10.1089/cmb.2019.0299, and the RECOMB-seq chapter 978-3-030-17083-7_14).
  * GitHub keeps the old objects reachable for a while; stripped data may remain
    downloadable via old SHAs until GitHub garbage-collects. Ask GitHub Support
    to run gc if that matters. This is public test data, so it is a tidiness
    issue rather than a disclosure one.

Recommended order (skip anything already done):

  1. Archive the data. ALREADY DONE if repo-slim-archive/ exists:
       tools/repo-slim/archive_data.sh
     Unlike isONcorrect's, this pulls the blobs out of HISTORY -- every file
     under test/ had already been deleted from HEAD, so .git was the only
     remaining copy. Verify the archive reads back before continuing:
       gzip -t  repo-slim-archive/raw/test/ccs.fastq.gz
       tar -tzf repo-slim-archive/raw/test/ENS_100k.fastq.tar.gz
       head -2  repo-slim-archive/raw/test/sample_alz_2k.fastq

  2. Back up the pre-rewrite history LOCALLY. Do not rely on a tag pushed to
     origin for this: filter-repo rewrites tags too, so the mirror push moves
     any such tag onto the rewritten commit. Worse, a tag left pointing at the
     old history would keep every stripped byte alive on the server AND get
     fetched by \`git clone\` (which pulls tags by default), undoing the whole
     exercise. Keep the old history off origin:
       git clone --mirror "$REPO_ROOT" ~/isONclust-preslim-backup.git

  3. Push the rewritten history. filter-repo removes the 'origin' remote from
     the mirror, so push to the URL rather than to a remote name:
       git -C "$MIRROR" push --force --mirror $ORIGIN

  4. Re-clone fresh and confirm:
       git clone $ORIGIN /tmp/isONclust-verify
       du -sh /tmp/isONclust-verify

  5. Put a fixture back. README.md and .travis.yml both reference one, so both
     are broken until this lands. Do NOT restore sample_alz_2k.fastq: it is
     blind to the minimizer window (PORTING.md, Finding 5). Use the replacement,
     which is 356 KB instead of 6.9 MB and discriminates --w:
       head -480 "$ISONCLUST_DATA/sirv/reads_10k_err7%.fastq" \\
         > /tmp/isONclust-verify/test/sirv_sim_120.fastq
       # expected sha256:
       #   85fd6f1d2d481365f07ffdd7bcdc305adab9a086102f8b59c84897f3eb38a35b
     ...and add a .gitignore for the junk this rewrite just removed, so it does
     not come back:
       printf '__pycache__/\\n*.py[cod]\\n.DS_Store\\n' >> /tmp/isONclust-verify/.gitignore

Step 3 is the irreversible one.
================================================================================
EOF

CLEANUP_WORKDIR=0
