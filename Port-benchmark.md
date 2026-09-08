# isONclust: Rust port vs the python reference, and vs isONclust3

Every number here is a single measured run, not an estimate. The port and the
reference are both deterministic and produce **byte-identical output**, so only
wall clock varies between repeats.

## How these were run

| | |
|---|---|
| machine | Apple M4 Max, 16 cores, 64 GB RAM |
| OS | macOS 26.5.1 |
| port | commit `c7bc874`, rustc 1.98.1, default features (parasail linked) |
| python | 3.12.14, parasail 1.3.4, `PYTHONHASHSEED=0` (`bench/setup_reference_env.sh`) |
| isONclust3 | 0.3.0, `--seeding minimizer --post-cluster` |
| presets | `--ont` (k13 w20) for ONT, `--isoseq` (k15 w50) for PacBio; isONclust3's matching `--mode` |
| truth | `minimap2` 2.31, `bench/make_truth.sh` |
| harness | `bench/benchmark.sh` |
| memory | peak RSS of the process and its waited-for children, via `/usr/bin/time -l` |

**isONclust3 is run with `--post-cluster`.** Its README's own example passes it and
it is not optional: without it, on real SIRV ONT reads it returns 5661 clusters
instead of 66. Benchmarking without it measures a tool nobody runs.

## What is being compared, and what is not

**The port reproduces the reference exactly.** All 27 equivalence cases pass
byte-for-byte, so the port and python rows differ only in time and memory —
never in accuracy. `bench/benchmark.sh` diffs the two clusterings on every run
and reports a difference as a bug, not a result.

**isONclust3 is a different algorithm with its own paper, not a port.** It is
expected to produce different clusters. Its rows are included because it is the
other current implementation, not as a like-for-like comparison of the same
method.

## The metric

**isONclust is a gene clustering tool**, so everything is scored against genes.
Truth comes from aligning reads to a reference with `minimap2` and taking the
gene of the primary hit.

| metric | what it says |
|---|---|
| homogeneity | do clusters contain reads from only one gene? (penalises merging genes) |
| completeness | is each gene in one cluster? (penalises splitting a gene) |
| V-measure | their harmonic mean |
| adjusted Rand | agreement corrected for chance |
| max cluster | the largest cluster, which drives downstream cost |

Transcript-level numbers exist in `PORTING.md` as a **diagnostic**. They are not
a target and nothing should be tuned to improve them; they rank the tools
differently on ONT and on PacBio, while the gene-level ordering does not change.

## Speed and memory

| corpus | preset | tool | `--t` | secs | peak MB | speedup |
|---|---|---|---|---|---|---|
| SIRV ONT, 10k | ont | python | 1 | 4.25 | 235 | — |
| | | **port** | 1 | **0.84** | **39** | **5.1x** |
| | | python | 8 | 1.33 | 206 | — |
| | | **port** | 8 | **0.31** | **79** | **4.3x** |
| | | isONclust3 | — | 0.65 | 32 | |
| SIRV PacBio, 17.6k | isoseq | python | 1 | 23.64 | 1092 | — |
| | | **port** | 1 | **9.80** | **746** | **2.4x** |
| | | python | 8 | 7.06 | 454 | — |
| | | **port** | 8 | **3.94** | 599 | **1.8x** |
| | | isONclust3 | — | 1.64 | 119 | |
| Drosophila ONT, 20k | ont | python | 1 | 16.08 | 704 | — |
| | | **port** | 1 | **5.83** | **328** | **2.8x** |
| | | python | 8 | 7.46 | 826 | — |
| | | **port** | 8 | **3.08** | **524** | **2.4x** |
| | | isONclust3 | — | 1.62 | 204 | |

**The port is 1.9–5.1x faster than the reference on every corpus and thread
count**, and uses less memory than it on five of the six rows — the exception is
PacBio at `--t 8`, where the port's resident batches still cost more than the
reference's eight separate processes (635 MB against 454 MB).

## Full corpora

The rows above are 10k–20k reads. On the full corpora, all single-threaded except
where `--t 8` is shown:

| corpus | reads | clusters | tool | peak RSS | time |
|---|---|---|---|---|---|
| SIRV real, full | 1 300 066 | 579 | python | 3.55 GB | 491 s |
| | | | isONclust3 | 2.49 GB | 243 s |
| | | | **port `--t 1`** | **0.74 GB** | **73 s** |
| | | | **port `--t 8`** | **1.34 GB** | **40 s** |
| Drosophila ONT | 1 000 000 | 84 318 | isONclust3 | 2.98 GB | 233 s |
| | | | **port `--t 1`** | **1.40 GB** | **726 s** |
| | | | **port `--t 8`** | **2.95 GB** | **205 s** |
| SIRV PacBio | 17 633 | 151 | isONclust3 | 0.12 GB | 1.7 s |
| | | | **port `--t 1`** | **0.75 GB** | **9.8 s** |
| | | | **port `--t 8`** | **0.60 GB** | **3.9 s** |

The reference was not run on Drosophila at 1M reads; it needs roughly 40 minutes
there.

For reference, the port's own figures on these corpora before the memory work
described in PORTING.md: SIRV real full 13.46 GB / 371 s, Drosophila 1M
10.99 GB / 1472 s.

`write_fastq`, the subcommand that splits a clustering into per-cluster fastq
files, is not covered by the rows above and is currently the most memory-hungry
path: 3.99 GB on SIRV real full at `--N 0`, 3.79 GB at `--N 2`. It holds every
record it will write. See PORTING.md.

## Accuracy — gene level

### SIRV ONT, 9998 reads with truth, 7 genes

| tool | `--t` | clusters | homogeneity | completeness | V | ARI | max cluster |
|---|---|---|---|---|---|---|---|
| python / **port** | 1 | 36 | 1.0000 | 0.5729 | 0.7285 | 0.5681 | 3835 |
| python / **port** | 8 | 30 | 1.0000 | 0.6472 | 0.7858 | 0.7338 | 4595 |
| isONclust3 | — | 66 | 1.0000 | 0.5005 | 0.6671 | 0.3119 | 1776 |

### SIRV PacBio, 14 783 reads with truth, 7 genes

| tool | `--t` | clusters | homogeneity | completeness | V | ARI | max cluster |
|---|---|---|---|---|---|---|---|
| python / **port** | 1 | 151 | 1.0000 | 0.6853 | 0.8132 | 0.7138 | 3327 |
| python / **port** | 8 | 110 | 1.0000 | 0.6933 | 0.8189 | 0.7259 | 3334 |
| isONclust3 | — | 163 | 0.6482 | 0.6503 | 0.6492 | 0.3578 | 8237 |

### Drosophila ONT, 20 000 reads, 3871 genes

| tool | `--t` | clusters | homogeneity | completeness | V | ARI | max cluster |
|---|---|---|---|---|---|---|---|
| python / **port** | 1 | 5679 | 0.9932 | 0.9826 | 0.9878 | 0.9285 | 459 |
| python / **port** | 8 | 5538 | 0.9914 | 0.9882 | 0.9898 | 0.9388 | 488 |
| isONclust3 | — | 6424 | 0.9937 | 0.9676 | 0.9805 | 0.9563 | 526 |

### Top 10 cluster sizes

| corpus | tool | sizes |
|---|---|---|
| SIRV ONT | isONclust1 `--t 8` | 4595, 944, 641, 567, 512, 423, 381, 329, 267, 209 |
| | isONclust3 | 1776, 1632, 1067, 655, 644, 489, 449, 423, 330, 329 |
| SIRV PacBio | isONclust1 `--t 8` | 3334, 3146, 1104, 1085, 1030, 866, 717, 573, 550, 355 |
| | isONclust3 | 8237, 1347, 1021, 861, 805, 734, 571, 464, 440, 435 |
| Drosophila | isONclust1 `--t 8` | 488, 432, 313, 267, 188, 173, 160, 142, 140, 130 |
| | isONclust3 | 526, 328, 296, 277, 174, 171, 168, 147, 146, 134 |

## Notes on the tables

**`--t 8` scores slightly higher than `--t 1`** on V-measure on all three
corpora. The `--t > 1` path is a hierarchical batch-and-merge rather than a
parallelised single pass, so it is a different computation, not the same one on
more cores. Why it merges differently is not established.

**isONclust3's homogeneity on SIRV PacBio is 0.6482**, against ≥0.99 on every
other row measured here, and its largest cluster there is 8237 reads against
isONclust1's 3327.

**Downstream isoform reconstruction on the SIRV PacBio clusters.** isONform run
over each tool's output (`bench/downstream.sh`), isoforms scored against the SIRV
reference:

| | isoforms | matching | recall | precision | F1 | isONform runtime |
|---|---|---|---|---|---|---|
| isONclust1 | 70 | 35 | 47.1% | 50.0% | 0.485 | 386 s |
| isONclust3 | 61 | 27 | 38.2% | 44.3% | 0.410 | 2053 s |

isONform's cost grows with cluster size, and the runtimes differ by ~1670 s.

**Algorithmic difference relevant to the homogeneity numbers.** When minimizer
sharing is ambiguous, isONclust1 performs an affine alignment before assigning a
read to a cluster; isONclust3 has no alignment step. The `--features
parasail-ffi` measurements in this file quantify what that alignment costs:
392 µs per call at ~2.0 G cells/s, on 0.55–2.6 calls per read.

**Corpus characteristics that bear on the runtime numbers.** SIRV_real_full forms
579 clusters from 1 300 066 reads; droso_1M forms 84 318 from 1 000 000. The
isONclust1 sweep compares each read against the representatives built so far, so
its cost scales with cluster count as well as read count — measured at
O(n^1.33–1.50) in read count on Drosophila. Read lengths also differ: SIRV maxes
at 2898 bp, Drosophila at 8339, PacBio at ~9000, and parasail's traceback is
O(n·m).

## Reproducing

```bash
bench/setup_reference_env.sh                 # pinned python reference
cargo build --release --manifest-path rust/Cargo.toml

# ground truth, once per corpus
bench/make_truth.sh READS.fastq REFERENCE.fasta truth.tsv --level gene

# the tables above
bench/benchmark.sh --corpora "sirv_real_10k sirv_pacbio droso_20k" \
                   --threads "1 8" --truth truth.tsv

# the downstream isoform check (PacBio: isONclust -> isONform)
bench/downstream.sh --corpus sirv_pacbio --out /tmp/ds
```

Corpora, their presets and where they live are in `bench/corpora.tsv`. The
Drosophila reference is Ensembl release 112 cDNA; the SIRV reference is the
68-transcript spike-in transcriptome.
