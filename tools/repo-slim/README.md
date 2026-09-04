# tools/repo-slim — shrinking the repository

The repo was ~492 MB to clone. **98.8% of that was committed test data** under `test/`. All source,
scripts and docs together are ~6 MB.

Measured on this history:

| | |
| --- | --- |
| total distinct blob bytes in history | 0.526 GB |
| stripped by these tools | 0.520 GB (98.8%) |
| retained source, scripts, docs | 0.006 GB |
| **rewritten `.git`** | **~1 MB** |
| **shrink factor** | **492×** |

## The three steps

Run them in this order. Each stops before doing anything irreversible.

```bash
tools/repo-slim/archive_data.sh     # 0. get the data out of git, into an archive
tools/repo-slim/analyze.sh          # 1. compute + review the removal list
tools/repo-slim/slim.sh             # 2. rewrite history in a scratch clone, verify
```

Only after all three pass do you force-push, and `slim.sh` prints the exact command. **Nothing here
pushes or uploads on its own.**

| Script | Does | Never does |
| --- | --- | --- |
| `archive_data.sh` | Extracts every `test/` blob from **history**, reassembles the split `ccs.fastq.gz`, writes a checksummed manifest | Upload, unless you pass `--upload` |
| `analyze.sh` | Writes `removal-paths.txt` and `analysis.txt` | Modify the repository at all |
| `slim.sh` | Clones a scratch mirror, rewrites it, verifies it | Touch your working repo, or push |

`removal-paths.txt` is generated, reviewable, and the single input to the rewrite. Read it first. It
came to 14 paths.

## Three things that made isONclust different from isONcorrect

These tools came from isONcorrect. Each of these needed the scripts changed, not just re-run.

**1. The data was already gone from `HEAD`.** Every file under `test/` had been deleted in the five
commits before this work began, so the isONcorrect version's premise — archive the working tree —
did not hold. `.git` was the only remaining copy of 519 MB. `archive_data.sh` now walks
`git log --all -- <path>` newest-first for the last commit that still had each file, and extracts the
blob. Run it *before* `slim.sh` or the data is gone for good.

**2. `test/ccs.fastq.gz` was committed as four `split` parts.** `part-aa` through `part-ad`, 365 MB
together — four blobs, not one file. Archiving them as-is would have preserved 365 MB of unusable
fragments, so `archive_data.sh` concatenates them and verifies the result with `gzip -t` before the
originals are destroyed. It passes; the file is 500 000 reads.

**3. There is a `0.0.4` tag whose tree contained stripped paths.** `filter-repo` rewrites tags, which
is what you want — but two failure modes are silent. The tag can vanish, or it can survive pointing
at an *unrewritten* commit, which would keep every stripped byte reachable on the server **and** get
pulled by `git clone`, which fetches tags by default, undoing the whole exercise. `slim.sh` now
asserts the tag set is preserved and that each tag's tree is clean.

## Requirements

`git-filter-repo` — `pipx install git-filter-repo` or `pip install git-filter-repo`. Point at an
unusual install with `GIT_FILTER_REPO=/path/to/git-filter-repo`.

macOS ships bash 3.2, which has no `mapfile` and errors on empty-array expansion under `set -u`.
These scripts are written for it. Do not add `mapfile`.

## The deduplication trap, and why it does not bite here

`git rev-list --objects` emits each blob exactly *once*, with a single path, so any file whose content
is byte-identical to another is invisible in that listing. Building the removal list from it silently
missed a file in isONcorrect: the rewrite ran, reported success, and left an 18 MB file in place.

Checked on this history: only three blobs are reachable at more than one path
(`modules/.DS_Store` == `scripts/.DS_Store`, both stripped;
`cemetary/cluster_parallel.py` == `modules/cluster_parallel.py`;
`modules/compute_shared_minimizers_probabilities.py` ==
`scripts/compute_shared_minimizer_probabilities.py`, both retained). So the trap does not bite —
but the removal list is still built from a tree walk (`git log --all --name-only`), because that is
the only way to *know* that.

The same trap applies to `git rev-parse HEAD:<missing-path>`, which prints the unresolved string to
stdout *and* exits non-zero. Existence checks use `git cat-file -e`.

## What the rewrite costs

- **Every commit SHA changes.** All 9 forks diverge permanently and cannot be fast-forwarded; anyone
  holding a clone (80 stargazers) must re-clone.
- **Commit-pinned links break**, including any in the two papers (10.1089/cmb.2019.0299 and the
  RECOMB-seq chapter 978-3-030-17083-7_14).
- **GitHub keeps old objects reachable for a while.** Stripped data may stay downloadable via old
  SHAs until GitHub garbage-collects; ask GitHub Support to run `gc` if it matters. This is public
  test data, so it is a tidiness question rather than a disclosure one.

Mitigation: back the old history up **locally**, not as a tag on origin. A tag pointing at old history
keeps every stripped byte alive on the server and gets fetched by `git clone`.

## Verification

`slim.sh` refuses to declare success unless all of these hold, and all of them did:

- commits touching source preserved exactly (80; the total legitimately falls 100 → 89 as data-only
  commits become empty and are pruned)
- 19 named source files have byte-identical blobs before and after
- known-stripped paths are absent from `HEAD`
- **no** path under `test/` or `modules/__pycache__/`, and no `.DS_Store`, survives anywhere in history
- the tag set is preserved and every tag's tree is clean

## Afterwards

Three things go back on top of the pushed history:

1. `test/sample_alz_2k.fastq` from `repo-slim-archive/raw/test/` — README links it as the install
   check and `.travis.yml` runs it, so both are broken until it lands. 6.9 MB; git stores content
   once, so re-adding costs exactly what keeping it would have.
   `sha256 35ed7203054c32484ee99a2ad26c29a12bf50758225d757dab53181a505d8751`
2. A `.gitignore` for the junk the rewrite removed, so it cannot come back.
3. `PORTING.md` and `bench/`, which is where the port starts.
