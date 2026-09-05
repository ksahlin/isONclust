# isONclust — Rust rewrite

## Goal

Port isONclust from Python to Rust. **Identical CLI, byte-identical output, faster, lower memory.**

The Python implementation in `isONclust` and `modules/` is the **normative reference**: when Rust and
Python disagree, Python is right until a human decides otherwise. Unlike the isONcorrect port, whose
specification was *accuracy*, this port's specification is **byte-identity**, and there are **no
deliberate divergences** until the port is exact. Improvements go in *Deferred improvements* and land
after exactness, each in its own commit.

The reference has one non-determinism defect of its own (*Finding 1*). **Decision taken: the Python
is left alone and the reference is pinned to Python ≥3.12**, where CPython's compensated `sum()`
closes the defect for us. The port therefore targets 3.12 semantics — exactly rounded summation —
and the goldens are valid only for that interpreter, which `bench/golden/manifest.tsv` records. The
one-line `math.fsum` fix that would make older interpreters agree is written up in *Deferred
improvements* and not applied.

### Why not just use isONclust3

The first question any reader has, and the README asks it for them — it opens by recommending
[isONclust3](https://github.com/aljpetri/isONclust3).

isONclust3 is already Rust, and it is faster and usually more accurate. It is **not a port**: it is a
different algorithm with its own paper, and it produces **different clusters**. That is the whole
distinction:

| | isONclust v1 (this port) | isONclust3 |
| --- | --- | --- |
| algorithm | greedy, quality-value based, `sorted(score)` single pass | different algorithm, own paper |
| output on the same input | **unchanged** | **changed** |
| what it offers | the same results, faster | better results |

So a byte-identical port of v1 is a genuinely separate deliverable. It makes the *existing* pipeline
faster without invalidating anything downstream of it — published results, snakemake workflows,
isONcorrect/isONform runs keyed to v1 cluster ids, and anyone who has already tuned thresholds
against v1 behaviour. isONclust3 asks all of those to be revalidated; this port asks for nothing.

Once the port is exact, isONclust3 becomes the comparison target rather than the alternative: **the
port and isONclust3 get measured against each other on accuracy and speed**, on the same corpora,
with the port standing in for v1. That comparison is only meaningful if the port is known-exact
first, which is why exactness comes before any optimisation.

The entry point must keep its exact current name, flags, and defaults:

- `isONclust` — clusters a fastq (or a `ccs.bam`/`flnc.bam` pair), writing `final_clusters.tsv`,
  `final_cluster_origins.tsv`, `sorted.fastq` and `logfile.txt` into `--outfolder`
- `isONclust write_fastq` — a subcommand that splits a clustering into per-cluster fastq files

## Branches

**Port work happens on `develop`, and `develop` merges to `master` when the port is exact.**

`develop` already existed and was not what it looked like. Measured before adopting it:

| | |
| --- | --- |
| `origin/develop` before this port | 80 commits, **9 behind `master`** |
| its only unique commit | `4822c01` "integrated medaka to polish consensus" |
| what `master` has that it lacked | `--k [4-9]` support (PR #13), fixes for issues #15 and #16, README updates |

So the pre-existing `develop` was a stale feature branch holding **exactly the work being dropped**,
and adopting it as-is would have silently reverted `--k [4-9]` and two issue fixes. It is therefore
re-created from `master`, which discards `4822c01`. That is intentional and is the same decision as
*Scope*'s: the consensus feature and its medaka polishing step go away.

`4822c01` is the only content lost, it is reachable in `~/isONclust-preslim-backup.git`, and it is
also the reason *Finding 2* needed a correction — medaka is genuinely implemented there, just never
released.

## Layout

| Path | Role |
| --- | --- |
| `isONclust` | Reference: CLI, the output writers, the `--consensus` driver, `write_fastq` |
| `modules/cluster.py` | Reference: minimizers, hit collection, the mapping and alignment decisions, `reads_to_clusters` |
| `modules/get_sorted_fastq_for_cluster.py` | Reference: quality scoring, filtering, the score sort, `sorted.fastq` |
| `modules/parallelize.py` | Reference: batching and the hierarchical merge driven by `--t` |
| `modules/consensus.py` | Reference: spoa subprocess, parasail alignment, reverse-complement detection. **Out of scope — the consensus feature is dropped.** Not a reference for anything the port implements |
| `modules/help_functions.py` | Reference: `readfq`, `cigar_to_seq`, `mkdir_p` |
| `modules/p_minimizers_shared.py` | Reference: a 2.5 MB Python literal — 59 628 rows of empirical minimizer-sharing probabilities. Data, not code |
| `rust/` | The port. Does not exist yet |
| `bench/` | The equivalence harness. `equivalence.sh`, `cases.tsv`, `diffsummary.py`, `setup_reference_env.sh`, `golden/` |
| `tools/repo-slim/` | The staged history-rewrite tooling. Already run; see *Repo hygiene* |
| `test/sirv_sim_120.fastq` | The committed smoke fixture. 356 KB, 120 simulated SIRV reads at 7% error, ground truth in every header. Linked from README, run by `.travis.yml` |
| `bench/corpora.tsv` | The corpus registry. Only the fixture is committed; the real corpora live under `$ISONCLUST_DATA` |
| `scripts/` | Paper experiment scripts. **Not part of the port. Do not modify.** `compute_cluster_quality.py` is the paper's accuracy scorer and will be reused for evaluation |
| `cemetary/` | Dead code the author kept. Not part of the port, not a reference for anything |

## Port status

Nothing is ported yet. This document is the reconnaissance, and it is the deliverable of the session
that wrote it.

| Stage | State | Verification |
| --- | --- | --- |
| reference environment | **done** | `bench/setup_reference_env.sh` builds it; parasail 1.3.4 from source on arm64 |
| reference runs on real data | **done** | droso 20k reads → 9255 clusters (1667 non-trivial) in 7.5 s |
| CLI contract captured | **done** | 14 cases in `bench/golden/cli/`, exit codes and stderr verbatim |
| output goldens recorded | **done** | 27 cases in `bench/golden/manifest.tsv` |
| determinism gate | **done, and it fails on Python ≤3.11** | `equivalence.sh seeds`; see *Finding 1* |
| interpreter decision | **taken: pin ≥3.12, reference unmodified for *Finding 1*** | recorded in the golden manifest |
| reference bug fixed upstream | **done — *Finding 9*'s `IndexError`** | its own commit on `master`; verified a no-op on 26 of 27 cases |
| corpora | **done** | 8 registered, 1 committed; the old fixture was replaced after being measured blind (*Finding 5*) |
| case matrix swept on all three corpora | **done** | `droso_20k` has 22 of 24 distinct results and **zero** unintended collisions — develop against it |
| repository slimmed and pushed | **done** | 492 MB → fresh clone at 4.8 MB |
| CLI parity | **done** | 23 differential cases green on the first run, 17 unit tests, clippy `-D warnings` and `cargo fmt` clean |
| `readfq` | **done** | 6 unit tests, incl. the no-trailing-newline path (*Finding 11*) |
| quality scoring + score sort (`sorted.fastq`) | **done** | `sorted.fastq` and `logfile.txt` byte-identical on 24 cases × 3 corpora, and on 10 configurations spanning ~257 000 reads. 12–15x faster |
| Python `str(float)` | **done** | `pyfloat.rs`; 39 997 values differentially checked against CPython |
| phred tables | **done, frozen** | *Finding 12*; 256 entries checked against the reference |
| `get_kmer_minimizers` | **done** | ~23 million minimizers identical across three corpora and six `(k, w)` settings, including *Finding 4*'s 1687 empty and 152 sub-k minimizers. `bench/dump_reference.py` + `equivalence.sh stage minimizers` |
| `get_all_hits` | **done, unit tests only** | not yet differentially verified — the replay oracle is fed the reference's own hit lists. `reads_to_clusters` is what will exercise it |
| `get_best_cluster` (the mapping decision) | **done** | ~95 000 recorded decisions identical across three real corpora, 54 000 of them assigning a cluster. Replayed from the live driver |
| Python `round(x, 2)` | **done** | `pyround.rs`; 100 602 values checked against CPython |
| empirical probability table | **done, frozen** | 437 KB blob generated from the 2.5 MB Python literal; *Finding 13* |
| `parasail_block_alignment` | **done** | 18 633 alignments identical — CIGAR *and* ratio — across four corpora and three of the four gap penalties. `parasail.rs` carried across from isONform |
| `get_best_cluster_block_align` | not started | the candidate walk around the alignment |
| `reads_to_clusters` (the driver) | not started | |
| output writers, cluster ordering | not started | |
| `parallelize.parallel_clustering` (`--t > 1`) | not started | semantic, not a speed knob; *Finding 3* |
| `write_fastq` | CLI only | parsing and its error paths are verified; the file-writing is not ported |
| `--consensus` (spoa + RC detection) | **dropped — will not be implemented** | *Scope*. `equivalence.sh dropped` asserts non-zero exit and the flag named |
| `--ccs`/`--flnc` (BAM input) | not started | needs a BAM reader; *Scope* |

## The pipeline, in one pass

One fastq in, one clustering out. There are no batches over reads in the v1 sense — the whole file is
sorted once and swept once.

1. **Score and filter every read** (`get_sorted_fastq_for_cluster`). For each read:
   homopolymer-compress the sequence and **drop it if `len(seq) < 2*k` or
   `len(hpol_compressed) < k`**. Compute `error_rate` as the mean per-base error probability from the
   quality string, and **drop the read if `10 * -log10(error_rate) <= --q`** (default 7.0). Then
   compute the score: `expected_number_of_erroneous_kmers` walks the quality string with a rolling
   product of per-base no-error probabilities, and
   `score = (1 - expected_errors/(len-k+1)) * (len-k+1)` — an estimate of how many error-free k-mers
   the read contains.

   Two different phred tables are in play and they are not interchangeable: `D` caps the per-base
   error probability at `0.79433`, `D_no_min` does not. The rolling product uses the **capped** one;
   `error_rate` uses the **uncapped** one.

2. **Sort descending by score** and write `sorted.fastq`. The score is **appended to every read
   accession** as `acc + "_" + str(score)`, and every later stage parses it back out with
   `float(acc.split("_")[-1])`. So `sorted.fastq` is simultaneously an output file and the input to
   everything downstream, and the float formatting is part of the contract. The sort is by score
   only, so ties fall back to input order (Python's sort is stable).

3. **Load the empirical probability table** (`p_minimizers_shared`). Keep rows where `k == --k` and
   `abs(w - --w) <= 2`, keyed by rounded error-rate pairs. This is the probability that two reads with
   given error rates share a minimizer, and it is what turns a gap between minimizer hits into a
   decision about whether the region is still "mapped".

4. **Initialise every read as its own cluster and its own representative.** Clustering only ever
   *merges* reads into an existing representative's cluster; it never creates a new cluster id.

5. **Sweep the reads in sorted order** (`reads_to_clusters`), highest score first. Per read:

   a. **Homopolymer-compress and take minimizers** (`get_kmer_minimizers`). Minimizers are chosen by
      **lexicographic order on the k-mer string** over a window of `w - k + 1` k-mers. Ties inside a
      window go to the **first** occurrence (`list(window).index`) — note this is the **opposite** of
      isONcorrect, which uses `rindex` and takes the last. Reads whose compressed length is under `k`
      are skipped here too, with a message on stdout.

   b. **Compute the compressed-read error rate.** Each homopolymer run contributes its single best
      quality character (`min` by error probability), and the mean error probability over that
      compressed quality string becomes the read's `error_rate` — a 7th element appended to its
      representative tuple. Note the representative tuple is 6 elements before this and 7 after, and
      the code branches on `len(...) == 7` to decide whether the work is already done.

   c. **Collect hits** (`get_all_hits`). For each minimizer, look it up in `minimizer_database`
      (k-mer string → set of representative ids) and, per hit representative, record the minimizer's
      index among the read's minimizers and its position in the compressed read. The read's own id is
      then deleted from the results.

   d. **Try to map** (`get_best_cluster`). Rank candidate representatives by
      `(number of shared minimizers, sum of hit positions, representative accession)` descending. If
      the best has fewer than `--min_shared` shared minimizers, give up immediately. Otherwise walk
      the ranking, stopping at the first candidate with fewer than `--min_fraction` of the top hit
      count. For each candidate, compute the probability that a *stretch between two consecutive
      minimizer hits* contains no shared minimizer purely by chance — a left-fold product of the
      per-minimizer error probability — and count the span between hits as **mapped** only when that
      probability is at least `--min_prob_no_hits`. The first candidate whose mapped fraction is
      **strictly greater than** `--mapped_threshold` wins.

   e. **If mapping failed but there were at least `--min_shared` shared minimizers, try aligning**
      (`get_best_cluster_block_align`). Only the candidates **tied at the top hit count** are
      considered. Sum the two reads' expected error rates, bin that into a gap-opening penalty
      (`5 / 4 / 3 / 2` at thresholds `0.01 / 0.04 / 0.1`), align with parasail
      `sg_trace_scan_16` (`match 2, mismatch -2, gap_ext 1`, falling back to `_32` on saturation),
      expand the CIGAR to a gapped alignment, then slide a window of `k` alignment columns and mark
      each window "aligned" if it contains at least `floor((1 - error_rate_sum) * k)` matches. The
      first candidate whose aligned fraction is at least `--aligned_threshold` wins.

   f. **Merge or become a representative.** If either step found a cluster, record
      `read → that cluster`. Otherwise the read stays its own representative and **all of its
      minimizers are added to `minimizer_database`**. This is what makes the sweep order-dependent:
      the database only ever grows, and a read is compared only against representatives that already
      exist, so processing order decides the outcome. That is why the score sort must be exact.

6. **Reassign.** Every read recorded in step 5f is moved into its target representative's cluster and
   its own now-empty cluster and representative entries are deleted. Merge targets are always
   representatives and representatives are never merge sources, so this is one level deep — there are
   no chains to follow. Worth pinning as an invariant in the port.

7. **Write output, ordered by `(cluster size, representative score)` descending.** Cluster ids in the
   output are *not* the internal ids: a fresh `output_cl_id` counts up from 0 in that sorted order.
   Ties in `(size, score)` fall back to dict order, which after step 6 is ascending internal id — so
   the port's ordering key is `(-size, -score, internal_id)`. Within a cluster, reads are written
   sorted by their score descending. The score suffix is **stripped** from accessions on the way out
   (`"_".join(acc.split("_")[:-1])`).

8. **Optionally build consensus sequences** (`--consensus`) — *described for completeness; this step
   is **out of scope** and the port will not implement it (see *Scope*).* Clusters at least
   `--abundance_ratio * total_reads` in size are handed to the `spoa` binary; the resulting centers
   are compared pairwise against each other's reverse complements with parasail, and any pair over
   `--rc_identity_threshold` identity is merged into the earlier one. Writes
   `consensus_references.fasta`.

### With `--t > 1` this is a different algorithm

`--t` does not parallelise the sweep above; it replaces it (`parallelize.parallel_clustering`).

Reads are cut into `--t` batches (by `--batch_type`: `total_nt` by default, or `nr_reads`, or
`read_lengths_squared`). Each batch is clustered **independently and in its own process**, producing
its own `minimizer_database`. The surviving representatives are then pooled, re-sorted by score,
regrouped by a threshold walk over their originating batch index, and clustered again — reusing the
minimizer database of the lowest-numbered batch in each group, and skipping the reads that already
built it. This repeats until one batch remains.

Consequences the port has to live with:

- **`--t` changes the answer.** On the smoke corpus `--t 1` gives 40 clusters and `--t 8` gives 35.
  On the retired `sample_alz_2k` the two agreed on 1235 of 1240 clusters and disagreed on 5 — so the
  size of the effect varies with the data, but it is never zero. Every `--t` value is its own
  equivalence case. A port that "parallelises with rayon" and expects `--t 1` output is not a port.
- **Parallel mode writes extra files.** One directory per merge iteration, each holding
  `pre_clusters.csv` and `cluster_origins.csv`. The iteration count depends on `--t`: 1 extra
  directory at `--t 2`, 2 at `--t 4`, 3 at `--t 8`. These are observable output and the harness
  diffs them.
- Processes are spawned with `mp.set_start_method('spawn')`, so workers do not inherit parent state.

## Scope: what gets ported

### Ported — inside the equivalence contract

`--fastq`, `--outfolder`, `--version`, `-h`/`--help`, `--k`, `--w`, `--q`, `--t`, `--d`,
`--ont`, `--isoseq`, `--min_shared`, `--mapped_threshold`, `--aligned_threshold`, `--min_fraction`,
`--min_prob_no_hits`, `--batch_type`, `--use_old_sorted_file`, and the `write_fastq` subcommand with
`--clusters`, `--fastq`, `--outfolder`, `--N`.

`--ont` is exactly `--k 13 --w 20`; `--isoseq` is exactly `--k 15 --w 50`. Both are pinned by their
own equivalence case, because a preset silently resolving to the wrong numbers is invisible in output
that agrees for other reasons.

### Deferred, not dropped

| Flag | Why deferred |
| --- | --- |
| `--ccs` / `--flnc` | BAM input via pysam. Needs a BAM reader (`noodles-bam`) and the flnc-within-ccs substring search that recovers quality values. Isolated from everything else — it only produces `read_array`. No corpus for it here; the archived `ccs.fastq.gz` is a fastq, not the BAM pair this path wants. |

### Dropped — the consensus feature, entirely

**Decision taken: `--consensus`, `--abundance_ratio`, `--rc_identity_threshold` and `--medaka` are
out of scope and will not be implemented.** They go away in the merge to main.

This is a deliberate narrowing of the tool, not a porting shortcut, and it is worth being clear about
what it buys and costs.

| | |
| --- | --- |
| **What it removes from the port** | the whole of `modules/consensus.py`: the spoa subprocess, a second parasail call site with different defaults (`opening_penalty=3` rather than the clustering path's error-rate-binned 2–5), the reverse-complement detection, and `consensus_references.fasta` as an output |
| **What it removes from the build** | the **spoa dependency, entirely**. The port needs no POA at all — no `spoars`, no carry-across from isONcorrect's `poa.rs`, and the harness stops caring whether `spoa` is on `PATH` |
| **What it removes from the contract** | one output file, and three flags' worth of equivalence cases |
| **What it costs** | anyone using `isONclust --consensus` loses it. It is a clustering tool; consensus of a cluster is isONcorrect's and isONform's job, and both do it better |

`--medaka` was already unusable on the release lineage — see *Finding 2*.

The port must **exit non-zero and name the flag** for all four. A pipeline that passes `--consensus`
and gets a zero exit with no `consensus_references.fasta` is worse off than one that fails loudly. A
generic "unrecognised argument" is not good enough either; the message should say the flag is not
supported. `bench/equivalence.sh dropped` asserts exactly that: **non-zero exit, and the flag named
in the output.**

Note this cannot be a recorded golden, because the reference does *not* do it — the reference
happily runs `--consensus` given spoa. It is the port's own contract and is asserted directly.

`--d 0` is not a flag to drop but an input to reproduce: it raises `ZeroDivisionError` and exits 1,
because `i % args.print_output` is evaluated before the truthiness guard that was meant to protect it.
It has a golden.

### Diagnostic only

`--d` (`print_output`) controls a progress table on **stderr** and does not change any output file.
Port it best-effort; matching the text is not part of the contract. `--verbose` does not exist here.

The reference also prints a good deal to **stdout** unconditionally — timings, the whole
`p_emp_probs` dict, cluster counts. Timings can never match, so stdout is **not** in the byte-identity
contract; the CLI goldens scrub timings and compare the rest.

## Determinism rules

The default code path is deterministic **on Python ≥3.12 only**, and that is not a property of the
algorithm. See *Finding 1* before reading the rest of this list.

- **`sum()` over floats is compensated from CPython 3.12 and a naive left-fold before it.** Four sites
  sum over `set(qual)`, whose iteration order is `PYTHONHASHSEED`-dependent. On ≥3.12 the interpreter
  makes the result order-independent; on ≤3.11 it does not. The port must reproduce the **exactly
  rounded** sum (Neumaier or equivalent), which then frees it to iterate in any order.
- **`reduce(mul, [p]*n, 1)` is not `p.powi(n)`.** It is a left-fold of `f64` multiplications starting
  from the integer `1`, and the port must fold the same way. This is `prob_all_errors_since_last_hit`,
  which decides what counts as mapped.
- **Minimizer ties go to the FIRST position in the window**, via `list(window_kmers).index(curr_min)`.
  isONcorrect takes the last. Do not carry that habit across.
- **`get_kmer_minimizers` can emit the same `(minimizer, position)` twice.** When the k-mer leaving
  the window merely *equals* the current minimizer by value, the minimum is recomputed and appended
  again, which can land on the same position. `get_all_hits` counts per minimizer index, so
  duplicates inflate hit counts. Reproduce it.
- **The sort in `get_best_cluster` is a total order**, because its third key is the representative's
  accession and accessions are unique. This matters: `minimizer_database` values are Python **sets of
  ints**, so the insertion order of `hit_clusters_hit_positions` is set-iteration order — but the sort
  key leaves no ties for that order to break, so it never reaches output. **Checked, not assumed**; it
  is the difference between needing a CPython-set model in Rust (isONform needed one) and not.
- **`sorted(reads)` is not used here.** Read order is the score-descending order written into
  `sorted.fastq`, and it is re-derived by parsing the score back out of each accession.
- **Cluster output ids depend on dict order.** `sorted(..., reverse=True)` is stable, so ties in
  `(size, score)` resolve to insertion order, which is ascending internal cluster id. Use an
  order-preserving map, or sort by `(-size, -score, id)`.
- **`write_fastq` iterates a `defaultdict` in insertion order**, which is the order cluster ids first
  appear in `final_clusters.tsv`. It decides nothing but the order files are created in.
- Dict iteration in the reference is insertion-ordered (Python 3.7+). Where iteration order feeds
  output, use an order-preserving map in Rust, not `HashMap`.

## Verification

Byte-identity is the acceptance criterion, and it is checked, not assumed.

```bash
bench/setup_reference_env.sh          # build the pinned reference env
bench/equivalence.sh env              # is it usable? is sum() compensated?
bench/equivalence.sh seeds            # does the reference agree with itself?  <-- run first
bench/equivalence.sh cli record        # capture the CLI contract
bench/equivalence.sh record           # record output goldens
bench/equivalence.sh verify           # run the port, diff against the goldens
```

**What counts as a difference.** Every file the tool writes, byte for byte:
`final_clusters.tsv`, `final_cluster_origins.tsv`, `sorted.fastq`, `logfile.txt`, and in parallel mode
`<n>/pre_clusters.csv` and `<n>/cluster_origins.csv` for each merge iteration. A file the port fails
to write and a file it writes that the reference does not are both failures.

`sorted.fastq` is **not** an intermediate to be skipped: the score is formatted into every accession
and parsed back out downstream, so it is an output and an input at once.

`logfile.txt` is included because it is the only place the error-rate distribution is observable, and
`error_rate` is precisely where the reference's determinism defect surfaces.

**The goldens are a manifest of hashes, not the files.** Recorded verbatim the 27 cases come to
318 MB, which has no business in a repository this exercise just took from 492 MB to 1 MB. Hashes are
enough to *fail* correctly; they cannot say *what* moved, so on a mismatch `verify` re-runs the
reference for that one case (about 2 seconds) and diffs properly. `bench/golden/manifest.tsv` is
137 KB and records the corpus hash, the interpreter version, and whether `sum()` was compensated —
because the goldens are only valid for the environment that produced them.

`bench/diffsummary.py` exists because a line diff is useless on these files:
`final_cluster_origins.tsv` carries the full read sequence and quality string in columns 3 and 4, so
one wrong float in column 6 prints four kilobytes. It reports which **column** moved, how many lines,
and the relative magnitude — and says outright when a float difference is below `1e-12`, i.e.
summation order rather than logic.

### The case matrix

27 cases in `bench/cases.tsv`, sweeping everything that can change output: `--k`/`--w` across four
settings plus both presets, `--q`, `--min_shared`, `--mapped_threshold`, `--aligned_threshold`,
`--min_fraction`, `--min_prob_no_hits`, `--t` at 1/2/4/8, all three `--batch_type` values, and
`write_fastq` at `--N` 0/2/10.

14 more cases in `bench/golden/cli/` pin the CLI contract: exit codes, the reference's own messages
verbatim, and the argparse behaviours a `clap` port will not reproduce by accident.

Several of those exit codes are wrong in an interesting way and are contract regardless:

| Case | Exit | Note |
| --- | --- | --- |
| no arguments | **0** | prints help and `sys.exit()` with no argument |
| `--ont --isoseq` together | **0** | prints "Arguments mutually exclusive" and exits 0 — a pipeline cannot detect this |
| `--flnc` without `--ccs` | **0** | same shape |
| `--fastq` with `--ccs` | **0** | same shape |
| `--w` < `--k`, or `--w` > 100 | 1 | one shared message for both |
| unknown flag | 2 | argparse's own |
| `--d 0` | 1 | `ZeroDivisionError` |
| `--medaka` | 1 | `UnboundLocalError` |

Also pinned: **argparse accepts any unambiguous prefix**, so `--outfold` works and `clap` will reject
it. And clap rewrites `field_name` to `--field-name`, so every multi-word flag needs an explicit
`long = "..."` — `--min_shared`, `--mapped_threshold`, `--aligned_threshold`, `--batch_type`,
`--min_fraction`, `--min_prob_no_hits`, `--use_old_sorted_file`, `--abundance_ratio`,
`--rc_identity_threshold`. The short-looking `--t`, `--d`, `--q`, `--k`, `--w` are **double-dash
single-letter** options, which is not what `clap`'s `short` produces either.

### The corpora

`bench/corpora.tsv` is the registry; `CORPUS` accepts a path or a name from it. Only the fixture is
committed.

| corpus | reads | truth | what it is for |
| --- | --- | --- | --- |
| `smoke` | 120 | per read (`@read_N_from_SIRV612`) | committed, 356 KB. CI, goldens, the README install check |
| `sirv_sim_err0` | 10 000 | per read | simulated, error-free. The high-quality extreme, where an integral score can appear and reach *Finding 6* |
| `sirv_sim_err7` | 10 000 | per read | simulated, 7% error. Paired with `err0`, this is the error-rate axis with everything else fixed |
| `sirv_real_10k` | 10 000 | transcriptome | real ONT. Discriminates `--w` sharply: 35 clusters vs 763 |
| `sirv_real_100k` | 100 000 | transcriptome | depth, for the `--t` sweep and profiling |
| `sirv_pacbio` | 17 633 | transcriptome | real PacBio CCS — the `--isoseq` path on data it was designed for |
| `droso_20k` | 20 000 | genome only | real ONT at transcriptome scale. **The only corpus that reaches *Finding 4*** |
| `droso_100k` | 100 000 | genome only | depth |

Goldens are corpus-specific and `manifest.tsv` records the corpus sha256, so changing `CORPUS` means
re-recording.

#### The fixture was replaced after it was measured

`test/sample_alz_2k.fastq` — 2500 PacBio CCS reads — was the corpus, is linked from the README, and
was run by `.travis.yml`. It was dropped, and the reason is *Finding 5*: at `--k 15` it produced
**byte-identical output for `--w 15` and `--w 50`**, so a port with a broken minimizer window would
have passed every case.

| corpus | `--w 15` vs `--w 50` at `--k 15` | distinct results across the 24 `main` cases | reads in the *Finding 4* zone |
| --- | --- | --- | --- |
| `sample_alz_2k` (retired, 6.9 MB) | **byte-identical** | 8 of 24 | 0 of 2409 |
| `smoke` (356 KB) | 70 vs 76 clusters | 15 of 24 | 0 of 120 |
| `sirv_real_10k` | 35 vs 763 clusters | not yet swept | 1 of 10 000 |
| `droso_20k` | not yet swept | not yet swept | **937 of 19 972** |

The replacement is 20× smaller and discriminates nearly twice as well. It is still only a smoke test:
`--q 0` and `--q 15` collapse onto the default on it, because simulated reads at a fixed error rate
have no quality spread. Two cases collide by design and should stay that way — `ont` == `k13w20`
pins that the preset resolves to `--k 13 --w 20`, and `t8` == `t8_total_nt` pins that `total_nt` is
the default `--batch_type`.

#### The matrices were swept, and `droso_20k` is the one to develop against

The full 27-case matrix, recorded on all three corpora. The question asked of each: **do any two
cases produce the same `final_clusters.tsv`?** A case that duplicates another is not a test.

| corpus | main cases | distinct results | unintended collisions | wall clock |
| --- | --- | --- | --- | --- |
| `smoke` (120 reads) | 24 | 15 | 2 | 12 s |
| `sirv_real_10k` | 23 | 17 | 2 | 1 m 41 s |
| **`droso_20k`** | **24** | **22** | **0** | 6 m 25 s |

Two collisions are intended on every corpus and must stay: `ont` == `k13w20` pins that the preset
resolves to `--k 13 --w 20`, and `t8` == `t8_total_nt` pins that `total_nt` is the default
`--batch_type`.

**`droso_20k` has no unintended collisions at all.** Every one of the 22 remaining cases produces a
result no other case produces, so every swept parameter is observable in output. That makes it the
corpus to develop the port against. What the weaker corpora hide:

| collision | on | why it hides something |
| --- | --- | --- |
| `default` == `q0` == `q15` == `mapped0.95` == `min_fraction0.5` == `min_fraction1.0` == `min_prob0.5` | `smoke` | simulated reads at one fixed error rate have no quality spread, so `--q` does nothing, and the thresholds never bind |
| `default` == `mapped0.95` == `min_fraction0.5` == `min_fraction1.0` | `sirv_real_10k` | the mapped fraction is near 1 on this data, so tightening the threshold changes nothing |
| `aligned0.9` == `k20w100`, `mapped0.3` == `min_prob0.01` | `smoke`, `sirv_real_10k` | coincidence at small scale |

`sirv_real_10k` is not redundant despite being weaker: it is the only corpus that found *Finding 9*,
because it is the only one whose reads are low-quality enough for `--q 12` to filter everything out.
`droso_20k`'s `--q 15` case passes cleanly. **Keep both.** This is method rule "weight the corpus by
its statistical power" and "three corpora, because one lies", arriving in the same afternoon.

Recording cost, for planning: 6 m 25 s for droso, and its manifest is 1.0 MB because `wf_N0` writes
9255 files. Neither is committed — goldens are per-corpus and recorded on demand.

### Stage-level oracles are required, not optional

*Finding 5* is the argument: end-to-end goldens on this corpus cannot observe `get_kmer_minimizers`
at all. A port with a wrong window would pass all 27 cases. So the port needs what isONcorrect and
isONform both needed — `bench/dump_reference.py` wrapping the reference **without modifying it**,
writing each stage's inputs *and* outputs in a stable line format, replayed from Rust and diffed
directly. Stages worth dumping, in dependency order:

1. ~~`readfq`~~ — **done**, verified through `sorted.fastq` itself
2. ~~`expected_number_of_erroneous_kmers` and the score~~ — **done**, same route
3. ~~the filter decisions~~ — **done**, same route
4. ~~`get_kmer_minimizers` — the full `(minimizer, position)` list per read, per `(k, w)`~~ —
   **done**, and this is the one that needed a dump: `bench/dump_reference.py --stage minimizers`
   against the port's `ISONCLUST_STAGE=minimizers`, both consuming the *reference's* `sorted.fastq`
   so a difference is the minimizer selection and not the sort
5. the compressed quality string and `error_rate`, per read
6. `get_all_hits` — the three dicts, per read
7. `get_best_cluster` — the ranking, the per-candidate mapped fraction, the winner
8. `parasail_block_alignment` — CIGAR, both gapped strings, alignment ratio
9. `reads_to_clusters` — the `read → cluster` decision and the database size, per read

End-to-end equivalence tells you *that* the port is wrong, never *where*.

**The dump oracle earns its keep, demonstrated.** Two plausible ways to get `get_kmer_minimizers`
wrong were introduced deliberately and both were caught on `droso_20k`:

| deliberate error | caught? | line count |
| --- | --- | --- |
| ties to the **last** position in the window (isONcorrect's rule, via `rposition`) | yes, all six `(k, w)` settings | **unchanged** |
| clamp the window to the sequence, i.e. "fix" *Finding 4* | yes, all six | **unchanged** |

Both produce exactly the same number of minimizers, so a check that compared counts — or a spot
check of the first few reads — would have passed. The comparison has to be the full ordered list.

The stage report also prints how many empty and sub-k minimizers each case actually exercised, so a
pass on a corpus that never reaches *Finding 4* says so rather than looking like coverage: the smoke
corpus reports `empty: 0, sub-k: 0` at every setting, and `droso_20k` reports up to 1687 and 152.

## Findings in the reference

Everything here was measured in the pinned environment, not inferred from reading. Each finding names
the corpus it was measured on, because that turned out to matter more than expected — see *Finding 5*.
Where a claim is "latent", it means the mechanism is confirmed but the corpus does not reach it, which
is a statement about the corpus and not a reason to skip the behaviour.

### Finding 1 — the reference does not agree with itself, and whether it does depends on the interpreter

**Stop-the-line.** Four sites compute an expected error count as
`sum([qual.count(c) * D[c] for c in set(qual)])`. `set(qual)` is a set of one-character strings, so
its iteration order is `PYTHONHASHSEED`-dependent — measured, the order differs on every seed.

Whether that reaches the result depends entirely on the CPython version:

| interpreter | `sum()` over floats | 400 orderings of one real quality string give |
| --- | --- | --- |
| 3.11.16 | naive left-to-right fold | **10 distinct sums**, spread 1.4e-15 |
| 3.12.14 | Neumaier compensated (gh-100425) | **1** — equal to the exactly rounded result |

And that reaches output. On Python 3.11, `--isoseq --t 1` on this corpus:

| file | across 5 seeds |
| --- | --- |
| `final_clusters.tsv` | stable |
| `sorted.fastq` | stable |
| `final_cluster_origins.tsv` | **810 of 1240 lines differ**, all in the `error_rate` column |
| `logfile.txt` | **differs**, in the highest and median error rate |

At `--t 8`, 405 lines differ. Every difference is a last-digit float artefact (relative 1e-16), which
is why `final_clusters.tsv` survives on *this* corpus — but the same perturbation feeds the `--q`
filter comparison, the `round(e, 2)` probability-table lookup, the four gap-penalty bins, and
`math.floor((1 - error_rate_sum) * k)`. Any of those can flip, and then the clustering changes.

So: **isONclust has been non-deterministic run-to-run on every interpreter it was released for**
(setup.py claims 3.4–3.7, `.travis.yml` tested 3.4–3.6), and it became accidentally deterministic on
3.12 through an unrelated CPython optimisation.

**Decision taken: pin Python ≥3.12 and leave the reference alone.** The goldens are therefore valid
only on 3.12+, which `bench/golden/manifest.tsv` records in its header, and `equivalence.sh env`
reports which behaviour the interpreter has. The port implements **exactly rounded** summation, which
matches 3.12 and — because the result is then order-independent — frees it to iterate in any order.
That is a real simplification: it is the difference between needing a model of CPython set iteration
in Rust (isONform needed one) and not.

What this decision costs, stated plainly: **the released tool stays non-deterministic for anyone on
Python ≤3.11**, and the port's contract now depends on an unrelated CPython optimisation rather than
on anything in the algorithm. If a user reports irreproducible `final_cluster_origins.tsv`, this is
why, and the fix below is the answer.

**The fix, not applied:** `math.fsum` instead of `sum`. Measured — `math.fsum` returns
`198.20590160207811` on both 3.11 and 3.12 for the same terms, order-independently on both, and that
is **the same value 3.12's `sum()` already returns**. So the fix changes nothing on 3.12, makes ≤3.11
agree with 3.12, and is available back to Python 2.6. Four call sites:
`get_sorted_fastq_for_cluster.calc_score_new`, `.fastq_single_core`, `.isoseq`, and
`cluster.reads_to_clusters`; plus `cluster.get_best_cluster_block_align`, which sums twice.

This would be method point 8 — fix the reference upstream once measured. It is **deferred by
decision**, not by oversight; it is listed in *Deferred improvements* and it is a one-line change per
call site whenever it is wanted.

`bench/equivalence.sh seeds` is the gate. It passes on 3.12 and **fails on 3.11**, which is how we
know it is a gate and not decoration.

### Finding 2 — `--medaka` crashes, and `--d 0` divides by zero

Both measured, both exit 1.

`--medaka` prints `Currently not implemented`, then falls through to `print("Saving references in:",
f.name)` — but `f` is only assigned in the *other* branch. `UnboundLocalError`.

**Correction, found when the branch strategy was settled:** that is true of `master` — the release
lineage, and what `0.0.6.1` on PyPI is — but **not** of `origin/develop`, whose single unique commit
(`4822c01`, "integrated medaka to polish consensus") actually implements it and removes the offending
`print`. So "nobody can be running this successfully" is wrong as stated: nobody running a *release*
can, and someone on `develop` can. The conclusion is unchanged, because the feature is now out of
scope entirely (*Scope*), but the reason is narrower than first written.

`--d 0` raises `ZeroDivisionError` at `if i % args.print_output == 0`. The `if args.print_output:`
guard that was presumably meant to protect it only gates the header line, several statements earlier.

### Finding 3 — `--t` is a semantic parameter wearing a performance parameter's clothing

`--t > 1` does not parallelise the sweep, it replaces it with a hierarchical batch-and-merge (see
*The pipeline*). Measured on this corpus, `--t 1` vs `--t 8`: both give 1240 clusters and 154
non-trivial ones, **1235 clusters are identical as read sets, and 5 differ each way**.

It also changes the file set: parallel mode writes `<n>/pre_clusters.csv` and
`<n>/cluster_origins.csv` per merge iteration — 1 iteration at `--t 2`, 2 at `--t 4`, 3 at `--t 8`.

The port cannot treat threading as an implementation detail. Every `--t` value is its own equivalence
case, and `--batch_type` is a third axis on top.

Also in `parallelize.py`: the `KeyboardInterrupt` handler calls `sys.exit()` and the module never
imports `sys`, so Ctrl-C raises `NameError` instead of exiting. Cosmetic, but it is what the reference
does.

### Finding 4 — `get_kmer_minimizers` reads past the end of the sequence, and emits sub-k minimizers

The initial window is built as `deque([seq[i:i+k] for i in range(w - k + 1)])` with no check that the
sequence is long enough. When the homopolymer-compressed read is at least `k` (the only guard) but
shorter than `w`, the later slices run off the end and Python silently returns short strings — and
then the empty string, which is lexicographically smallest and therefore wins.

Measured on an 18 nt sequence:

| `k` | `w` | minimizers returned |
| --- | --- | --- |
| 15 | 50 | `[('', 18)]` — one **empty-string** minimizer |
| 15 | 20 | `[('ACGTACGTACGTAC', 4)]` — a **14-mer**, for `k=15` |
| 13 | 20 | `[('ACGTACGTACGTA', 0)]` — correct |

**This is live on real ONT data, at the default settings.** Measured on `droso_20k`
(`--k 15 --w 50`, which is the parser default and also `--isoseq`):

| | count |
| --- | --- |
| reads surviving the upstream filter | 19 972 |
| with compressed length in `[k, w)` — the overrun zone | 937 |
| emitting an **empty-string** minimizer | **445** |
| emitting a **sub-k** minimizer | 152 |

At `--ont` (`--k 13 --w 20`) the zone shrinks to 55 reads: 0 empty-string, 18 sub-k. So this is
substantially a consequence of running ONT data at the *default* `--w 50` rather than at `--ont`,
which is easy to do because `--w 50` is the default with no preset flag at all.

**And the default `--min_shared` masks it.** A read that produces an empty-string minimizer produces
**exactly one** minimizer, so it can never reach `--min_shared 5` and is forced to be a singleton.
Measured: all 445 of them are singletons, and the largest cluster containing one has size 1. The
`''` key accumulates all 445 ids in the database and is never useful.

Lower the threshold and the collision appears:

| | clusters | largest cluster |
| --- | --- | --- |
| `--min_shared 5` (default) | 9 255 | 320 |
| `--min_shared 1` | 4 424 | **449** |

449 ≈ the 445 empty-string reads merged into one cluster on the strength of a shared empty string.
So the bug is real, reachable through a documented flag, and hidden by a default — the most durable
kind. The cost at the default is milder but not nothing: 445 reads that might have clustered are
guaranteed not to.

Reproduce all of it. Do not fix it in the port; it is in *Deferred improvements*.

### Finding 5 — the original corpus could not see the minimizer window at all

The most important finding for the harness. It retired the corpus and it is why stage oracles are
mandatory.

On `sample_alz_2k` at `--k 15`, `--w 15` and `--w 50` produced **byte-identical `final_clusters.tsv`,
`final_cluster_origins.tsv`, `sorted.fastq` and `logfile.txt`**. Verified by hand outside the harness,
not just as a golden collision.

That is not because the parameter is inert:

- minimizer density differs 19× — 1543 minimizers at `w=15` vs 80 at `w=50`, on one 1557 nt compressed read
- the probability table differs in **all 225 entries** (`w=15` gives `p=0.809246` where `w=50` gives
  `0.790523` at `e=(0.01, 0.01)`)

The clustering was simply insensitive to both on data that easy — which meant **a port with a broken
minimizer window would have passed every one of the 27 end-to-end cases.** Only 8 of the 24 `main`
cases were mutually distinct.

**Resolved by replacing the corpus.** The same comparison elsewhere:

| corpus | `--w 15` | `--w 50` | distinct results, 24 `main` cases |
| --- | --- | --- | --- |
| `sample_alz_2k` (retired) | — | — | byte-identical; 8 of 24 |
| `smoke` (committed, 120 reads) | 70 clusters | 76 clusters | 15 of 24 |
| `sirv_real_10k` | **35 clusters** | **763 clusters** | not yet swept |

The replacement fixture is 20× smaller than what it replaced and separates nearly twice as many
cases. On real SIRV ONT the parameter is not subtle at all: 35 clusters against 763.

Two lessons, and the second is the one that generalises:

1. **A corpus being real is not the same as a corpus being informative.** `sample_alz_2k` was real
   data, and useless for this question.
2. **The blindness was invisible from inside the harness.** 27 green cases looked like 27 checks.
   What exposed it was asking a question the harness does not ask — "do any two of these cases
   actually differ?" — which took one loop over the manifest. Worth running whenever the matrix or the
   corpus changes; a case that duplicates another is not a test, and the two collisions that remain
   (`ont` == `k13w20`, `t8` == `t8_total_nt`) are kept deliberately because each pins a default.

### Finding 6 — the score is a Python float repr, and it round-trips through a filename-shaped string

`sorted.fastq` headers are `acc + "_" + str(score)`, and downstream code recovers the score with
`float(acc.split("_")[-1])`. So Python's float formatting is in the byte-identity contract, and Rust's
`{}` does not match it:

| value | Python | Rust `{}` |
| --- | --- | --- |
| `1234.0` | `1234.0` | `1234` |
| `1e-5` | `1e-05` | `0.00001` |
| `1e16` | `1e+16` | `10000000000000000` |

**Latent on this corpus:** of 2408 distinct scores, **zero** are integral and **zero** use exponent
notation — they are all ordinary values like `1000.065598053912`, where Python's shortest-round-trip
repr and Rust's agree. A score is `(1 - errors/(len-k+1)) * (len-k+1)`, which is integral exactly when
the read has no expected errors, so a high-quality read reaches it. The port needs Python-compatible
float formatting; it just will not find out from this corpus.

Note also that accessions are split on `_` and the reference rewrites spaces to `_` when reading
fastq headers, so a read whose *original* name ends in `_<number>` is ambiguous. The corpus's
accessions contain `_` throughout (`..._s1_p0/3382/2212_58_CCS_strand=-;...`), and it works only
because the score is appended last. Do not "improve" the delimiter.

### Finding 7 — two quality tables, differing by one cap

`get_sorted_fastq_for_cluster` defines both `D` (per-base error probability, capped at `0.79433`) and
`D_no_min` (uncapped), and uses them for different things: the rolling no-error product uses the
capped table, `error_rate` uses the uncapped one. `cluster.py` defines its own copy of the **capped**
one as `phred_char_to_p`. Three tables, two distinct value sets, and the choice is load-bearing
because `error_rate` gates the `--q` filter.

The cap bites only where the uncapped probability exceeds `0.79433`, which is every character below
ASCII 34 — in practice just `!` (phred 0, uncapped probability `1.0`, capped to `0.79433`). Phred 1
(`"`) computes to `0.7943282347242815`, just under the cap, so `min` leaves it alone. **This is live
on the corpus: 121 reads contain `!`.** So the two tables genuinely disagree on real input, and using
the wrong one shifts `error_rate` for exactly those reads.

### Finding 8 — `--use_old_sorted_file` truncates the logfile and skips writing it

`get_sorted_fastq_for_cluster.main` opens `logfile.txt` with mode `w` as its *first* action, then
returns early if the sorted file already exists and the flag is set. So the logfile is emptied and
never rewritten. Observable, and the harness will see it.

### Finding 9 — `--q` above the corpus's quality crashed with an `IndexError`. Fixed in the Python

`error_rates[int(len(error_rates)/2)]` on a sorted list takes the upper middle element for even
lengths, with no averaging. Reproduce it.

The more serious half was written up as "needs a fixture" and then **the `sirv_real_10k` sweep hit
it**: if every read is filtered out, `min_e = error_rates[0]` raises `IndexError` on an empty list.
Measured, on real SIRV ONT reads:

| `--q` | exit | reads passing |
| --- | --- | --- |
| 7 (default) | 0 | 9 972 |
| 8 | 0 | 8 407 |
| 9 | 0 | 4 433 |
| 10 | 0 | 369 |
| 11 | 0 | **2** |
| 12 | **1** | 0 |
| 15 | **1** | 0 |

So a user who raises `--q` to be stricter on ONT data gets `IndexError: list index out of range`
instead of "no reads passed the filter". `--q 12` is an entirely reasonable thing to ask for.

**And it is corpus-dependent, which is the point of having more than one corpus.** The same
`--q 15` case runs fine on `droso_20k` — those reads are higher quality — and would have been
recorded as a clean pass had `sirv_real_10k` not been swept. One corpus lies.

**Fixed in the reference**, on `master`, as its own commit — this is method point 8, and the first
time it has been applied in this port. The statistics are guarded; the tool now says
`Error: no reads passed the quality filter (--q 12.0).` on stderr, notes the same in `logfile.txt`,
and exits 1.

The fix was scoped and then verified rather than assumed:

| | |
| --- | --- |
| exit status | **unchanged at 1** — an unhandled exception already exited 1, so nothing checking the exit code sees a difference |
| `sorted.fastq` | **unchanged**, still created and still empty |
| `logfile.txt` | 0 bytes → one line of explanation. The only output byte that moved |
| smoke corpus, all 27 cases | **byte-identical** before and after |
| `sirv_real_10k`, all 27 cases | **only `q15` differs, and only in `logfile.txt`**; the other 26 byte-identical |

That last row is the check worth copying: a bug fix in the reference is a behaviour change, so it has
to be measured against the full matrix on a corpus that reaches it, not just against the case that
crashed. The port therefore targets the fixed behaviour, and the goldens were re-recorded — which on
the smoke corpus changed nothing at all.

### Finding 14 — the mapping oracle has two blind spots, and they are worth naming

`equivalence.sh stage mapping` replays ~95 000 recorded `get_best_cluster` calls and they all match.
Before trusting that, three deliberate errors were introduced. **One was caught and two were not:**

| deliberate error | caught? | why |
| --- | --- | --- |
| candidate ranking sorted ascending instead of descending | **yes**, both settings | changes which candidate wins |
| `reduce(mul, ...)` replaced by `p.powi(n)` | **no** | see below |
| `mapped_ratio >= threshold` instead of `>` | **no** | see below |

Neither miss is a broken oracle; both are the corpus being unable to distinguish the cases, and the
measurements say so precisely:

* **The fold and `powi` genuinely differ** — in the last bit, for **27 783 of the 32 400 `(p, n)`
  combinations this corpus actually reaches** (e.g. `p=0.17140336964776648, n=3` gives
  `0.005035679329970428` against `0.005035679329970429`). They still produce the same *decisions*,
  because the value is only ever compared against `min_prob_no_hits`, and one ULP flips that only
  when it sits exactly on the threshold. That never happened in 95 000 calls.
* **`>` versus `>=`** is observable only when `mapped_ratio` equals `mapped_threshold` exactly.
  Across 4098 distinct ratio values on `sirv_real_10k`, **none** was exactly 0.7.

Both are implemented faithfully anyway, because exactness is the specification, and both are pinned
by unit tests instead. The honest summary is that this oracle checks the *decision* strongly and the
*arithmetic* weakly — and a green run should not be read as more than that.

A third thing the corpus decides: **both simulated corpora map nothing.** `Passed mapping criteria`
is **0** for the smoke fixture and 0 for `sirv_sim_err7` — every assignment goes through the
alignment fallback — against 8144 for `sirv_real_10k`, 12 967 for `sirv_pacbio` and 3950 for
`droso_20k`. Running this stage on the committed fixture proves nothing at all, so the harness prints
how many reads were assigned and warns when the answer is zero.

### Finding 13 — fifteen CLI-valid `(k, w)` settings crash on an empty probability table

`p_emp_probs` is built by keeping table rows where `k == args.k` and `abs(w - args.w) <= 2`. The
table's `w` values step by 5 for each `k`, starting at `k` — so `k=6` has 6, 11, 16, … 96, and
nothing within ±2 of `w=99` or `w=100`. The dict comes out empty and the first lookup dies:

```
KeyError: (0.01, 0.01)
```

Measured: **15 `(k, w)` combinations that pass the CLI's own `k <= w <= 100` validation** reach it,
including `--k 6 --w 100`, `--k 7 --w 100` and `--k 16 --w 100`. All are at the top of the window
range, and all are reasonable things to ask for.

Cheap to fix properly — the check belongs beside the existing window validation, where it can say
which settings are available instead of failing several seconds into a run. In *Deferred
improvements*; the port reports it rather than raising a `KeyError`.

### Finding 11 — a fastq without a trailing newline crashes the tool

`readfq` truncates each line with `l[:-1]`, which removes the last character
unconditionally rather than stripping a newline. On the final line of a file that
does not end in one, the running quality length never reaches `len(seq)`, the loop hits EOF, and the
reference falls through to `yield name, (seq, None)`.

The caller does not guard it. Measured, on a two-read fastq with no trailing newline:

| | exit | result |
| --- | --- | --- |
| with trailing newline | 0 | 2 reads clustered |
| **without** | **1** | `TypeError: 'NoneType' object is not iterable` |

Files get truncated, hand-edited and generated without a final newline all the time, so this is easy
to hit and gives a traceback rather than an explanation.

**Fixed in the reference**, on `master`, as its own commit. The fix is minimal: a line that ends in
`\n` is chomped exactly as `l[:-1]` did — carriage returns included, so CRLF files are untouched —
and only the final unterminated line changes. The quality-length counter follows the same rule
instead of assuming `len(l) - 1`.

Verified a no-op on everything that already worked: **0 differing manifest rows** across the full
27-case matrix on the smoke corpus, on `sirv_real_10k` and on `droso_20k`. All three end in a
newline, which is both why they were unaffected and why the bug survived this long.

A quality string that is genuinely shorter than its sequence still yields `None` and still crashes
the reference. That path is untouched, and the port reports it rather than reproducing the
traceback.

Worth noting how this was found: the first version of the port's unit test asserted that the last
*quality value* was dropped, which is what reading the code suggests. Asking the reference showed the
whole quality string is dropped instead. The test was wrong, not the port.

### Finding 12 — Rust's `powf` and CPython's `**` differ by one ULP, and only real ONT reads notice

The sorting stage was byte-identical on the smoke corpus, on both simulated SIRV corpora, and on real
PacBio CCS — and **differed on every real ONT corpus**, in the 13th significant digit of the score:

```
@read_37_..._strand=+_652.3743164321687     reference
@read_37_..._strand=+_652.3743164321701     port
```

The cause is one entry of the 95-entry phred table. `10 ** (-(ord(c) - 33) / 10.0)` for `%`
(phred 4, i.e. `10 ** -0.4`):

| how it is computed | bits |
| --- | --- |
| CPython `10 ** -0.4` | `3fd97a967f7524b2` |
| Rust `10f64.powf(x)` at runtime | `3fd97a967f7524b3` |
| C `pow(10.0, -0.4)` against this machine's libm | `3fd97a967f7524b3` |
| Rust `powf` on a literal (LLVM constant-folds it) | `3fd97a967f7524b2` |

One value of 95, one ULP — and `%` is an ordinary ONT quality character that simulated reads and
PacBio CCS reads happen not to contain. The error then compounds through
`expected_number_of_erroneous_kmers`'s rolling product, which is why it surfaces at 1e-14 rather than
1e-16.

**The deeper finding is about the reference, not the port:** CPython disagrees with this machine's own
C library, so the reference's value depends on which libm its interpreter was built against. The
reference is not portable here either, and neither would a port be that computed the table.

Resolved by freezing both tables as generated constants (`rust/src/phred.rs`), exactly as the
argparse text is frozen. That makes the port reproducible on every platform and identical to the
pinned reference. `rust/tests/phred_oracle.rs` re-checks all 256 entries against the live reference,
so a changed environment fails a test instead of drifting silently.

**This is what the corpus sweep was for.** Four of seven corpora agreed while the port was wrong. Had
`sample_alz_2k` still been the only corpus — PacBio CCS, no `%` — this would have shipped.

### Finding 10 — harness bugs found while building the harness, both silent

Recorded because both are the shape of mistake this project will make again.

**A vestigial `shift` ate the first argument.** `cli_case` had a leftover
`local target="$1"; shift`, so `cli_case version --version` invoked the reference with *no*
arguments. It recorded exit 0 and a 4558-byte stdout, which is exactly what `--version` failing to
`--help` looks like if you do not check the content. Two goldens were wrong and the harness reported
`14 passed`.

**`cases.tsv` was space-aligned instead of tab-separated.** With `IFS=$'\t'`, `read -r name entry
args` put the entire line into `$name`; `$args` was empty; **every one of the 27 cases ran with no
arguments and recorded the same default output.** The harness reported 27 passes, and the file count
per case was a plausible 10. The tell was the case names in the log having the arguments concatenated
onto them — visible, and easy to read past.

**A case that writes no files was silently skipped.** `wf_N10` asks for clusters of 10+ reads and
the smoke corpus has none, so the case legitimately produces zero files — and therefore contributed
zero rows to the hash manifest. `verify` looked up its expected exit code, found nothing, reported
`no golden for wf_N10` and moved on: **neither a pass nor a failure**. 26 of 27 cases failed and the
27th was invisible. The port could have done anything there. Fixed by writing one `__meta__` row per
case unconditionally, carrying the exit code and the file count, so "no golden" now means genuinely
absent and is itself a failure.

**Two goldens contained run-varying data, so they could never match anything.** Found by noticing
that `git status` was dirty immediately after a re-record. `bench/golden/cli/medaka/stdout` held a
`tempfile.mkdtemp()` path (`/var/folders/.../tmpybxe1c8f`), and `manifest.tsv` held its own
`# recorded:` wall-clock line. A golden like that is not a check — it is a permanent failure, and a
permanent failure trains you to ignore the harness. The temp path is now scrubbed to `<TMPDIR>`
alongside the timings, and the manifest carries no timestamp (git records when it was committed; what
matters for validity is the corpus hash and the interpreter, which are still there).

All three are now guarded:

| guard | catches |
| --- | --- |
| `check_cases` | a case line that is not tab-separated into three fields |
| `equivalence.sh stable` | a golden containing a timestamp, temp path, PID or duration — records twice and diffs |
| `equivalence.sh verify` | demonstrated on **a single digit** changed in one field of one line out of 1240, reporting file, line, column and magnitude |
| `equivalence.sh dropped` | demonstrated against three stand-ins: refuses-and-names passes, accepts-and-ignores fails, refuses-vaguely fails |
| `equivalence.sh seeds` | demonstrated to fail on Python 3.11 |

**A failing case aborted the whole run, twice.** `cli_case` ended with a diff whose non-zero status
became the function's return value, and under `set -e` that killed the script: the first CLI verify
reported one failure out of 28 cases and the other 27 had never run. The same shape came back in the
dump-based stage helpers, where a trailing `[[ ... ]] && info ...` left status 1 whenever the
condition was false — so a *passing* mapping run printed its results and then exited before the
`==> N passed` summary, with a non-zero status. Both fixed with an explicit `return 0`.

This is the single most recurrent mistake in the harness, and it fails in the most misleading
direction available: the run looks like it stopped because something went wrong.

Four harness bugs, every one of them reporting something other than the truth — three reported
passes, the fourth reported a single failure and hid 27 unrun cases. The gates below exist because of
them, and each was verified to fail before being trusted.

The general rule this earns is in *Method* below: a harness that has never failed has not been tested,
it has only been run.

## The reference environment is not `pip install -r requirements.txt`

`bench/setup_reference_env.sh` builds it, and the reason it is a script rather than a line of
documentation is that parasail fights back on Apple Silicon.

parasail publishes **no macOS arm64 wheel** — 1.3.4 ships `macosx_10_9_x86_64` only — and it is not
on conda-forge or bioconda under any name (`parasail`, `parasail-python`, `libparasail`: all absent).
So it must build from source, and its `setup.py`:

1. **Prepends `/usr/bin` to `PATH`** before probing for build tools, so macOS's `m4 1.4.6` always
   beats a newer `m4` on your `PATH`. Having decided `m4` is too old, it downloads `m4-1.4.17` from
   `ftp.gnu.org`, builds it, and then fails to run it. The reported error is
   `RuntimeError: autoreconf -fi failed`, which names neither `m4` nor the download. Setting `$M4`
   is the way out, because `autoreconf` honours it.
2. **Requires the tools spelled `glibtoolize` and `glibtool`** on Darwin — the Homebrew names.
   conda-forge's libtool installs them as `libtoolize` and `libtool`, so the probe fails on a
   perfectly good libtool and says nothing at all about which tool was missing.

Two symlinks and one environment variable. Pinned versions in
`bench/env/resolved-*.txt`; the built environment is python 3.12.14, pysam 0.24.0, parasail 1.3.4.

**`spoa` is not needed at all.** It was the one remaining external tool, required only by
`--consensus`, which is now out of scope (*Scope*). `equivalence.sh env` no longer looks for it.
That leaves the reference environment at exactly two dependencies — parasail and pysam — and pysam
only for the deferred BAM path.

**Build both interpreters.** `bench/setup_reference_env.sh 3.11 isonclust-ref-311` is what makes
*Finding 1* reproducible rather than a claim, and it is the only way to check that the determinism
gate still fails when it should.

## Repo hygiene

Done and **pushed**, before any port work existed, because it ends in a force-push and a port branch
in flight would have made that worse.

`tools/repo-slim/{archive_data.sh,analyze.sh,slim.sh}` came across from isONcorrect and were adapted.
The measured result: **492 MB → 1 MB, 492× smaller**, 14 paths stripped, 98.8% of all blob bytes.

| | |
| --- | --- |
| total distinct blob bytes in history | 0.526 GB |
| stripped | 0.520 GB (98.8%) |
| retained source, scripts, docs | 0.006 GB |
| rewritten `.git` | ~1 MB |

What was stripped: `test/ccs.fastq.gz.part-{aa,ab,ac,ad}` (365 MB), `test/old_sorted_ens_100k.fastq.tar.gz`
(74 MB), `test/ENS_100k.fastq.tar.gz` (73 MB), `test/sample_alz_2k.fastq` (6.9 MB),
`test/isonclust1.out`, four committed `.pyc` files under `modules/__pycache__/` (1.15 MB, wrong for
any interpreter but 3.6), and two `.DS_Store` files.

**Three ways isONclust's situation differed from isONcorrect's**, all of which needed the tooling
changed rather than just re-run:

1. **The data was already gone from `HEAD`.** Every file under `test/` had been deleted in the five
   commits before this work started, so `archive_data.sh`'s premise — archive the working tree —
   did not hold: `.git` was the only remaining copy. It now extracts blobs from **history**, walking
   `git log --all -- <path>` newest-first for the last commit that still had each file.
2. **`test/ccs.fastq.gz` was committed as four `split` parts.** Four blobs, not one file. Archiving
   them individually would have preserved 365 MB of unusable fragments, so `archive_data.sh`
   reassembles and verifies with `gzip -t` before the originals are destroyed. It passes; the file is
   500 000 reads.
3. **There is a `0.0.4` tag, and its tree contained stripped paths.** `filter-repo` rewrites tags, so
   this works — but two silent failure modes exist: the tag can vanish, or it can survive pointing at
   an *unrewritten* commit, which would keep every stripped byte reachable on the server **and** get
   pulled by `git clone`, which fetches tags by default. `slim.sh` now asserts the tag set is
   preserved and that each tag's tree is clean. Both pass.

`slim.sh` verified 19 source files byte-identical before and after, that commits touching source are
preserved (80), that no stripped path survives anywhere in history, and the tag checks above — all
green — and then stopped, because pushing is a separate human decision. Authorised and executed:

```
+ 6119457...4822c01 develop -> develop (forced update)
+ 004e74f...795de23 master  -> master  (forced update)
+ d55cabc...10685cd 0.0.4   -> 0.0.4   (forced update)
```

Confirmed by a fresh clone:

| | before | after |
| --- | --- | --- |
| `git clone` total | ~500 MB | **4.8 MB** |
| `.git` | 492 MB | **2.2 MB** |
| commits | 100 | 89 |
| tag `0.0.4` | present | present, tree clean |
| stripped paths anywhere in history | 14 | **0** |
| the seven reference source files | — | byte-identical |

The deduplication trap that made isONcorrect's first removal list wrong does not bite here — checked,
only three blobs in all of history are reachable at more than one path — but the removal list is still
built from a tree walk, because that is the only way to *know* that.

### The data is not deleted from GitHub, and cannot be

The intent was to remove it permanently rather than archive it publicly. That is not achievable, and
it is worth being exact about what was and was not accomplished. Measured immediately after the push:

| | result |
| --- | --- |
| `git fetch origin <old master sha>` | **refused** — `upload-pack: not our ref` |
| `GET /repos/ksahlin/isONclust/commits/004e74f` | **still returns the commit** |
| `GET /repos/.../git/blobs/1546863d…` (the 6.9 MB fixture) | **still returns 6 865 440 bytes** |
| `GET /repos/.../contents/test/…?ref=<old sha>` | 404 |
| GitHub's reported repository size | **still 501 MB, unchanged** |
| `test/` on the default branch of `sukses24/isONclust`, `nextgenusfs/isONclust` | **all four `ccs.fastq.gz` parts plus `ENS_100k.fastq.tar.gz`, live** |

So: normal git access no longer serves the old history, but the objects are still there and still
addressable by SHA through the API. Two reasons, and only one of them is fixable:

1. **GitHub has not garbage-collected.** It runs on its own schedule. A support request asking them
   to `gc` the repository would drop the now-unreachable objects — worth opening if this matters.
2. **The forks legitimately still contain the data on their live default branches**, and GitHub
   shares object storage across a fork network. Those are other people's repositories. Nothing done
   here, and nothing GitHub Support can do short of detaching the repository from its network, makes
   that data unreachable — and detaching breaks the fork relationships.

This is public test data behind a published paper — simulated CCS reads and a public PacBio sample —
so it is a tidiness question rather than a disclosure one. **What was actually achieved is the thing
that mattered: the repository clones in 4.8 MB instead of 500 MB.** Full removal was never available.

**Archive disposition: local only.** `repo-slim-archive/` (gitignored) holds all 495 MB, verified
readable, with a checksummed manifest. Not published to Zenodo or a Release. It is the only
convenient copy, so do not delete it without deciding that the data is genuinely unwanted — though
per the table above, the forks are an inconvenient backstop.

After the push, four things go back on top: `test/sirv_sim_120.fastq` (356 KB — README links it and
`.travis.yml` runs it, so both are broken until it lands), a `.gitignore` for the junk the rewrite
removed, `bench/` and `tools/repo-slim/`, and this file. **The old fixture was deliberately not
restored** — see *Finding 5*.

## Deferred improvements

Nothing here may land before the port is byte-identical, and each needs its own commit.

### Known bugs in the reference

Ordered by how much they matter.

1. **`sum()` over a `set` — non-determinism.** *Finding 1*. `math.fsum` at five call sites.
   **Deferred by decision**, not by oversight: the reference is pinned to Python ≥3.12, where the
   interpreter closes the defect. Applying it would make ≤3.11 agree with 3.12 at no cost on 3.12
   (measured), so it is available whenever someone on an older interpreter reports irreproducible
   `final_cluster_origins.tsv`. If it is applied, re-record the goldens.
2. **`get_kmer_minimizers` reads past the end of the sequence.** *Finding 4*, and **live on real
   data**: 445 of 19 972 `droso_20k` reads emit an empty-string minimizer at the default
   `--k 15 --w 50`. The guard checks `len(hpol) < k` but the window needs `w`. Correct fix is to
   clamp the window, or to skip reads shorter than `w`; either changes clustering for short reads, so
   it is a behaviour change with an accuracy question attached. `droso_20k` is the corpus to measure
   it on, and `--min_shared 1` is the setting that makes the damage visible (largest cluster 320 →
   449).
3. **`--medaka` crashes.** *Finding 2*. Either implement it or remove the flag.
4. **`--d 0` divides by zero.** *Finding 2*. Move the modulo inside the truthiness guard.
5. ~~**`error_rates[0]` on an empty list.**~~ **Fixed** on `master`, its own commit, and merged to
   `develop`. *Finding 9* has the measurements. The port targets the fixed behaviour.
6. **The "mutually exclusive" and "needs both" validation paths exit 0.** A wrapper script cannot
   detect them. Changing them to exit non-zero is right and is a breaking change for anyone who
   depends on the current behaviour, which is presumably nobody.
7. **`parallelize.py` calls `sys.exit()` without importing `sys`.** *Finding 3*. Ctrl-C during
   parallel clustering raises `NameError`.
8. **A file handle is opened per cluster in the `--consensus` loop and closed only inside the
   `if`.** Clusters below the abundance cutoff leak a handle and truncate `reads_tmp.fq`. Harmless to
   output; will bite on a large cluster count. **Moot for the port** — `--consensus` is dropped — and
   recorded only because it is still a live bug in the Python for as long as the Python ships it.
9. **`detect_reverse_complements` compares against centers already merged away**, and only the last
   element takes the `i == len - 1` branch. The merging is order-dependent in a way that is probably
   not intended. **Also moot for the port**, same reason.

### The aligner: is parasail the bottleneck, and would WFA2 do?

**Measured before answering.** Profiling `reads_to_clusters` on `droso_20k` (19 938 reads, 10 971
alignment calls), single core:

| | cumulative | share |
| --- | --- | --- |
| `reads_to_clusters` total | 26.63 s | 100% |
| `get_best_cluster_block_align` | 12.42 s | 47% |
| ├─ `parasail_block_alignment` | 10.79 s | 41% |
| │  └─ **`parasail.sg_trace_scan_16` itself** | **4.60 s** | **17%** |
| │  └─ `cigar_to_seq` | 1.19 s | 4% |
| `get_kmer_minimizers` | 3.30 s | 12% |
| `get_all_hits` | 2.15 s | 8% |
| `get_best_cluster` | 1.37 s | 5% |

So in the reference, **parasail is 17% and the Python wrapped around it is 23%** — the CIGAR
decoding, the gapped-string expansion and the rolling match window cost more than the alignment.
Replacing the aligner and keeping the wrapper would chase the smaller half.

That is a statement about the *Python*, and it inverts after porting: the wrapper cost largely
disappears in Rust while parasail's does not, so the aligner's share of what remains goes up. Method
point 5 — profile before optimising, **and re-profile after**. The number to act on is the one
measured against the finished port, not this one.

**Would WFA2 do?** isONform already answered most of this, and the answer is more interesting than
"it is faster":

* **WFA2's ends-free mode is not parasail's semi-global** (isONform, finding 40). parasail
  *maximises a score* with a positive match reward; WFA2 *minimises a penalty* with `match = 0`.
  With all four ends free the empty alignment costs nothing and is therefore always optimal — on two
  identical 20 bp sequences WFA2 returns `DDDD…IIII`, aligning nothing. Even with the ends bounded it
  pays nothing for terminal matches and declines them: 48 bases shaved off two identical 200 bp
  sequences.
* **isONform's reconciliation is the reusable part.** `src/wfa.rs` bounds the free ends, greedily
  extends the aligned core outward over matching pairs, then scores the resulting columns with
  **parasail's own rules**. WFA2 only chooses the columns; the `Scoring` prices them. That is what
  keeps the arithmetic checkable, and it is directly portable here because this port already carries
  isONform's `parasail.rs` and its `Scoring`.
* **It was a real improvement there**, not a trade: on droso, +14 FSM isoforms over both Python and
  the faithful port, and 4.7x faster (isONform, finding 55).

**On gap penalties and exon differences** — the specific worry. Both libraries use the same affine
model, `open + (L - 1) * ext`, so WFA2 is not inherently worse: configured with isONclust's
parameters it prices a gap identically. Neither is *well suited* to exon-scale gaps, though. With
`gap_ext = 1` a skipped 500 bp exon costs about 500, which at `match = 2` needs 250 matching bases
just to break even, so a genuine exon difference tends to look like a bad alignment rather than a
gap. What WFA2 offers that parasail does not is **two-piece affine** (a second, cheaper regime for
long gaps), which is exactly the model this problem wants. That is the interesting question here, and
it is not primarily a speed question at all: it would change which reads cluster together, probably
for the better, and on ONT data the alignment path decides **10 309 of 19 938 reads** on
`droso_20k` — the majority.

**And the exact port is 15x slower than the C library.** Measured on `droso_20k`: replaying the
reference's own 11 241 alignments takes **67 s** in the port against **4.6 s** for
`parasail.sg_trace_scan_16` (dump parsing is 0.02 s of that, so it is all alignment). `parasail.rs`
is an exact scalar dynamic program — it exists to be *right*, and reproducing parasail's tie-breaking
bit for bit is what it is for. It is not, and was never meant to be, fast.

That inverts the picture above and is the single most important performance fact so far: **on ONT
data the port will be slower than the reference until this is addressed**, because the alignment path
decides the majority of reads. The 12–15x won on the sorting stage does not pay for it.

isONform's answer was `simd.rs` — block-aligner behind a parasail-shaped call — but that is an
*approximation*: isONcorrect measured it as "optimal score on 1400/1400 recorded alignments; differs
only in which equally-optimal path it reports, changing ~0.8% of reads". Under a byte-identity goal
that is a divergence, not an optimisation, so it cannot be adopted while exactness is the
specification. It is the obvious candidate the moment that changes.

**Which corpora reach which gap penalties**, measured across the recorded alignments — the penalty is
chosen per comparison from the two reads' summed error rates, so this is really a statement about
read quality:

| corpus | alignments | opening penalties reached |
| --- | --- | --- |
| `smoke` | 101 | **5** only |
| `sirv_real_10k` | 1 806 | **2** only |
| `sirv_pacbio` | 5 485 | **4** only |
| `droso_20k` | 11 241 | 2 (3 823), **3** (7 247), 4 (171) |

No single corpus reaches all four. `droso_20k` covers three, and only the near-perfect reads in
`smoke` reach 5 — which is the one case simulated data is *good* for. **`smoke` + `droso_20k`
together cover all four**, and that pair is the minimum any aligner comparison must run on. Running
one corpus alone would leave a quarter of the parameter space unmeasured, and `sirv_real_10k` alone
would leave three quarters.

**So the work, in order:**

1. Finish the port and **re-profile**. The 17% above is a Python measurement and will not survive.
2. If the aligner is then the bottleneck, port isONform's `wfa.rs` alongside its `parasail.rs`,
   keeping parasail as the exact baseline. Gate the swap on **clustering verdicts**, not on
   alignment scores (isONform, finding 41).
3. Two-piece affine is a **separate** experiment from the WFA2 swap, and a behaviour change either
   way. Measure it alone, against an exact baseline, on `droso_20k` where the alignment path is
   load-bearing — isONform's finding 55 exists because stacking divergences inverted the sign of an
   earlier conclusion.

Every step here is a deliberate divergence and needs its own commit and its own note, per the goal:
byte-identity first, improvements after.

### Memory: the whole dataset is resident, and 2-bit encoding is the obvious win

**Yes, every read is in memory at once, and more than once.** Confirmed by reading and by measuring:

- `get_sorted_fastq_for_cluster` builds `read_array` holding `(acc, seq, qual, score)` for every
  surviving read, sorts it, and writes `sorted.fastq`.
- `isONclust.main` then **reads that file back** into a second `read_array` of
  `(i, b_i, acc, seq, qual, score)`.
- `single_clustering` copies each entry into `representatives`, and `reads_to_clusters` replaces
  each with a 7-tuple carrying the same `seq` and `qual` again.
- In parallel mode (`--t > 1`) the batches are **pickled to worker processes**, so the peak is
  multiplied by the number of cores.

Measured on `droso_100k` — 99 547 reads, a 136 MB fastq:

| | peak RSS |
| --- | --- |
| reference, full run at `--t 1` | **1476 MB** (~11x the input) |
| the port's sorting stage alone | 608 MB (~4.5x) |

**The proposal: pack nucleotides two bits each** (`A=00, C=01, G=10, T=11`), which cuts sequence
storage 4x. Worth doing, and worth knowing exactly what it costs before it lands:

- **It breaks the byte-identity contract on any dataset containing `N` or another non-ACGT
  character**, if those are mapped to a pseudo-random nucleotide. This is not hypothetical: `N` is
  observable in at least two places. `get_kmer_minimizers` picks minimizers by **lexicographic order
  on the k-mer string**, and `'N'` (ASCII 78) sorts between `'G'` (71) and `'T'` (84) — so an `N`
  changes which k-mer wins a window. And parasail's matrix is built for `"ACGT"`, scoring anything
  else as 0 rather than as a match, which changes the alignment path.
- **Every corpus in `bench/corpora.tsv` is pure ACGT** — checked, zero non-ACGT bases in the smoke
  fixture, `sirv_real_10k`, `droso_20k` and `sirv_pacbio` — so the change would be *measurably*
  lossless on everything currently tested. That is a reason to be careful rather than reassured: it
  means the test suite cannot see the divergence, which is *Finding 5*'s lesson again. A corpus with
  `N`s has to be built before this lands.
- **A third option avoids the contract break entirely:** three bits per base, or 2-bit plus a
  sparse side-table of exception positions. `N` is rare in ONT and PacBio output, so a side-table
  costs almost nothing and keeps the port exact. Prefer that unless measurement says the extra
  indirection is expensive.
- **Quality strings are the other half of the footprint and cannot be 2-bit** — they carry 40+
  distinct values. They also cannot simply be dropped after scoring: `reads_to_clusters` needs the
  quality string again to compute the homopolymer-compressed error rate, and
  `get_best_cluster_block_align` recomputes `poisson_mean` from the *full* quality string for every
  candidate it considers. Caching one float per read instead of retaining the string is likely the
  larger and safer win, and it is behaviour-neutral.

Order of work, when it comes: build an `N`-containing corpus; cache the per-read quality statistics;
then pack sequences, with the exception table; measure each separately.

### Performance and structure, once exact

- **The minimizer database is `k-mer string → set of ids`.** Every lookup hashes a `k`-length string.
  Interning k-mers to a 2-bit-packed `u64` (k ≤ 32 covers the full CLI range, 4–30) removes the string
  hashing from the innermost loop. Behaviour-neutral: the map is insertion-ordered but its iteration
  order is proven not to reach output (*Determinism rules*).
- **`get_all_hits` allocates three `defaultdict`s per read** and appends to two of them per hit.
  One pass over a reusable buffer keyed by candidate id is the obvious replacement.
- **`prob_all_errors_since_last_hit` builds a full list of `n+1` products per candidate**, where the
  per-minimizer probability is a single value repeated. It is a prefix-product of a constant, so the
  whole list is `p^gap` — but it must still be computed as a left-fold of that many multiplications to
  stay bit-identical, so the win is allocation, not arithmetic.
- **`qual.count(char_)` inside a comprehension over `set(qual)`** is O(len × distinct) per read. One
  histogram pass is O(len). This is the same code *Finding 1* fixes, so do both at once — and note the
  histogram makes the summation order canonical for free, which is the actual fix.
- **`--t` spawns processes and pickles the whole read array per batch.** In-process threads with a
  shared read array is the port's natural shape, and it must reproduce the batch boundaries exactly.
- **Profile before optimising, and re-profile after.** No profile has been taken yet. The reference
  runs this corpus in 1.8 s, which is too small to profile; the archived `ccs.fastq.gz`
  (500 000 reads) is the right size and is in `repo-slim-archive/`.

## Method

Carried over from the isONcorrect and isONform ports. The full versions, with the measurements behind
each point, are in those repositories' `PORTING.md`.

1. **CLI parity first**, locked by unit tests. Argument names, defaults, validation order, stderr text
   and exit codes. Watch that clap rewrites `field_name` to `--field-name` — every multi-word flag
   needs an explicit `long = "..."`, and this reference has nine of them plus five double-dash
   single-letter flags.
2. **Differential oracles, not end-to-end tests.** Wrap the reference without modifying it; dump each
   stage's inputs *and* outputs in a stable line format; replay pure functions from Rust and diff
   stateful ones. End-to-end equivalence tells you *that* something is wrong, never *where*.
   *Finding 5* makes this mandatory here rather than advisable.
3. **Dump from the live driver too, not only from a standalone dump binary.** Oracles replay recorded
   reference inputs, so a port following a different trajectory still passes them. The one real bug in
   the isONcorrect port passed every oracle and was caught only by diffing a dump taken from the
   running driver.
4. **Build a real corpus before trusting anything.** Simulated and spike-in data gave the wrong answer
   twice in isONcorrect. Here the corpus is real but weak (*Finding 5*), which is a different failure
   with the same remedy.
5. **Profile before optimising, and re-profile after.** Bottlenecks move. Instrument at stage
   granularity, and *remove* sub-stage instrumentation after reading it: on millions of calls the
   timers themselves distorted the table.
6. **Measure, do not reason, about performance.** Recorded null results from isONcorrect: caching an
   edit-distance pattern was slower, reusing the POA engine was worth nothing, 4-bit DP cells were
   slower than 8-bit, hoisting a hash lookup out of a loop was worth nothing.
7. **Set up CI on day one, on Linux *and* macOS, x86_64 *and* arm64.** Three defects in isONcorrect
   existed only because every local check ran on one machine. Note for this repo: `.travis.yml` is
   dead (Travis, python 3.4–3.6, and it runs the fixture that was just deleted from `HEAD`), so CI is
   new work, not a migration.
8. **Fix reference bugs upstream once measured, rather than reproducing them.** The largest accuracy
   win in the isONcorrect port was a one-character fix to the *Python*. Each such fix is its own
   commit, with goldens re-recorded. *Finding 1* is this port's first and it is a blocker.

### Four more rules, earned in the isONform port after isONcorrect's list was written

* **A conclusion at one depth is not a conclusion.** In isONform, `delta_len=3` fixed one corpus at
  20 000 reads and did nothing at 5 000 or 10 000; a median spoa backbone won at one depth and lost at
  another. Sweep the depth before writing the sentence.
* **Weight the corpus by its statistical power.** Six depths of a 68-transcript panel put every arm
  within ±2 transcripts, non-monotonically — noise. A 561-cluster annotated corpus separated them
  without ambiguity. Do not pick a default from a corpus that cannot distinguish the options. *This
  port's corpus cannot distinguish `--w 15` from `--w 50`* (*Finding 5*); that is the same rule
  arriving early.
* **A mechanism confirmed on one case is not a prediction about the population.** In isONform a
  mechanism genuinely caused one instance's over-extension, and turning it off still did not help the
  corpus.
* **Measure one divergence at a time, against an exact baseline.** A whole run of findings was
  measured with other divergences active, and one of them turned out to have the sign inverted. This
  port's "no deliberate divergences until exact" rule exists for the same reason.

### Two more from isONform, about the checks rather than the conclusions

* **A differential harness is only as good as the question it asks.** isONform's CI compared the port
  against the reference after the port's *default* had moved to optimised semantics — so every
  differential step was comparing two programs that are meant to differ and calling it a failure,
  while `cargo test` stayed green the whole time. Changing the default changed the question without
  changing the harness. Corollary earned here, in *Finding 10*: **a harness that has never failed has
  not been tested, it has only been run.** Both harness bugs found while building this one reported
  passes. Deliberately break the port and confirm the harness notices — `verify` has been shown to
  catch one digit in one field of 1240 lines, and `seeds` has been shown to fail on Python 3.11.
* **A check that runs on one machine measures that machine.** isONform's first CI run found two things
  no amount of local checking would have: a clippy lint that did not exist on the dev machine's
  toolchain but did on CI's floating `stable`, and a filename whose case only matters on Linux
  (`isonform_parallel` vs `isONform_parallel` — macOS is case-insensitive, so the wrong spelling
  resolved everywhere it had ever been run). A third thing hid alongside them: a job targeting a
  retired GitHub runner **queues indefinitely rather than failing**, so the matrix looked like four
  platforms and was three, silently. Relevant here immediately — this reference's entry point is
  `isONclust` with three capitals, and everything so far has been measured on one arm64 Mac.

## Working agreements

**Report results, not progress.** While anything is running, say nothing — no interim findings, no
partial tables, no "here is what I am about to do". One work stint gets one reply, in three parts and
no more:

1. **what was measured**, in a sentence or two;
2. **the numbers**, as a table;
3. **what to do next**, as a proposal.

Never restate a conclusion already given. Corrections are one sentence: what is now true, not a
retrospective on the error. Detail goes in this file, not in the reply — the reader will ask if they
want more.

- **Never run `git commit`, `git push`, or any history-changing git command.** Leave finished work in
  the tree, say what changed and why, propose the commit message, and let a human commit it. This held
  throughout the isONcorrect and isONform ports and holds here — and it holds especially for the
  history rewrite in `tools/repo-slim/`, which is a force-push over 9 forks and 80 stargazers.
- Don't "improve" the algorithm while porting. A behaviour change and the port must not land in the
  same commit; an intentional divergence needs its own commit and a note here. **This port's
  specification is byte-identity, so there should be no intentional divergences at all** until it is
  exact.
- **"Behaviour" means observable output, not internal representation.** Different containers, dropping
  provably-dead entries, arena allocation, 2-bit-packed k-mers — none of that is a behaviour change if
  the emitted bytes are identical. What is not free: iteration order where it reaches results,
  tie-breaking, arithmetic and rounding. In this reference specifically: minimizer tie direction,
  summation order, the left-fold of products, float formatting, and cluster-id assignment order.
- When you spot a possible improvement, write it into *Deferred improvements* and move on.

## First steps, in order

1. ~~Slim the repository before any port work exists.~~ **Done, verified, pushed.** 492 MB → a 4.8 MB
   fresh clone. *Repo hygiene*.
2. ~~Get the reference running in a pinned environment.~~ **Done**, on both 3.12 and 3.11.
   `bench/setup_reference_env.sh`.
3. ~~Check `PYTHONHASHSEED` sensitivity before recording anything.~~ **Done, and it found a defect.**
   *Finding 1*.
4. ~~Decide what to do about it.~~ **Done: pin Python ≥3.12, reference unmodified.** The `math.fsum`
   fix is written up and deferred. Recorded in the golden manifest so the constraint travels with
   the goldens.
5. ~~Capture the CLI contract and record output goldens.~~ **Done.** 14 CLI cases, 27 output cases.
6. ~~Replace the corpus.~~ **Done**, and it was necessary: the old fixture was blind to the minimizer
   window (*Finding 5*). `bench/corpora.tsv` registers 8; one is committed.
7. **Commit the work.** The push replaced origin's history, so everything in this session's tree is
   untracked relative to it. Proposed commits are listed below.
8. ~~Sweep the full case matrix on `droso_20k` and `sirv_real_10k`.~~ **Done, and it found a
   reference crash.** `droso_20k` discriminates every swept parameter (22 of 24 distinct, zero
   unintended collisions) and is the corpus to develop against; `sirv_real_10k` is weaker but is the
   only one that exposed *Finding 9*, where `--q 12` and above crashed with an `IndexError` on real
   ONT reads. Keep both. *Finding 9* is now fixed in the Python.
9. **Port the CLI**, locked by unit tests and the 14 differential cases. Nine multi-word flags need
   explicit `long = "..."`; five are double-dash single-letter; argparse prefix abbreviation is live.
10. ~~Work outward from the leaves: `readfq`, the quality scoring and the score sort.~~ **Done.**
    `sorted.fastq` and `logfile.txt` are byte-identical across 24 cases on three corpora and on ten
    configurations spanning ~257 000 reads, and the stage is 12–15x faster (droso_100k: 7.90 s →
    0.51 s). Verified with `bench/equivalence.sh stage sort` rather than a dump/replay oracle: this
    stage's output *is* a file the tool writes, so diffing it directly is simpler and stronger.
    **`get_kmer_minimizers` is next and does need `bench/dump_reference.py`**, because its output
    never reaches a file — and *Finding 5* is the argument that end-to-end goldens cannot see it.
11. **CI on Linux and macOS, x86_64 and arm64, on day one.** Method point 7. `.travis.yml` is dead
    (Travis, Python 3.4–3.6) and should be replaced, not migrated. Note the entry point is `isONclust`
    with three capitals — the exact trap that broke isONform's CI on its first ext4 filesystem.
12. Then `--t > 1`, then the BAM path. **Not `--consensus`** — it is dropped (*Scope*), which also
    means the port never needs spoa.
13. **Then, and only then, compare against isONclust3** on accuracy and speed, with the port standing
    in for v1. That comparison is the point of the exercise and it is only meaningful once the port is
    known-exact.

## Proposed commits

Nothing in this session was committed — see *Working agreements*. `.reconcile/RECONCILE.md` holds a
runnable script for all of it.

**First, reconcile the working copy.** The force-push replaced origin's history, so this clone's
`master` is disjoint from `origin/master` (99 commits vs 88, no shared SHA) and its `.git` still
holds 490 MB of unreachable objects. That step was deliberately not run: `git reset --hard` is a
history-changing command, and it would discard the two *tracked* files this session modified
(`README.md`, `.travis.yml`). Both are backed up as a patch, verified to apply cleanly.

**Then three commits on `master`, because the slimming left `master` broken.** `README.md` links
`test/sample_alz_2k.fastq` as the install check and `.travis.yml` runs it, and that file no longer
exists in the repository. This is repo hygiene, not port work, so it does not belong on a port
branch:

| # | Branch | Contents | Message |
| --- | --- | --- | --- |
| 1 | `master` | `.gitignore` | `Ignore build artifacts and OS cruft stripped from history` |
| 2 | `master` | `tools/repo-slim/` | `tools: add the staged history-rewrite tooling` |
| 3 | `master` | `test/sirv_sim_120.fastq`, `README.md`, `.travis.yml` | `test: replace the install fixture with one that discriminates --w` |

Commit 3 is the one carrying a claim: the README's stated cluster counts (35/21 at the default
`--t 8`, 40/24 at `--t 1`) are measured and reproducible. `.travis.yml` is updated only so it stops
naming a file that no longer exists — it is dead CI either way (Travis, Python 3.4–3.6).

**Then `develop`, re-created from `master`, carrying the port.** See *Branches* for why `develop`
has to be re-created rather than adopted:

| # | Branch | Contents | Message |
| --- | --- | --- | --- |
| 4 | `develop` | `bench/` | `bench: add the equivalence harness, corpora registry and goldens` |
| 5 | `develop` | `PORTING.md` | `Add the Rust port plan, reconnaissance and findings` |

Pushing `develop` is a **force-push that discards `4822c01`**, the medaka implementation. That is the
intent (*Scope*), but it is worth doing deliberately rather than as a side effect. Two ways:

```bash
# (a) keep the commit reachable by name, then discard the branch
git tag archive/develop-medaka origin/develop && git push origin archive/develop-medaka
git push --force origin develop

# (b) just discard it -- it stays in ~/isONclust-preslim-backup.git either way
git push --force origin develop
```

(a) costs one tag and a few kilobytes, and means nobody has to know about a local backup mirror to
find the medaka work again. Recommended.
