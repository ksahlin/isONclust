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

## What is being compared

Memory and runtime of the port is compared against the Python implementation.
Accuracy is identical since **the port reproduces the reference exactly.** 

Memory, runtime, and accuracy is compared against isONclust3 (which is a 
different algorithm, also in Rust).

## The metric

**isONclust is a gene clustering tool**, so everything is scored against genes.
Truth comes from aligning reads to a reference with `minimap2` and taking the
gene of the primary hit.

| metric | what it says |
|---|---|
| homogeneity | do clusters contain reads from only one gene? (penalizes merging genes) |
| completeness | is each gene in one cluster? (penalizes splitting a gene) |
| V-measure | their harmonic mean |
| adjusted Rand | agreement corrected for chance |
| max cluster | the largest cluster, which drives downstream cost |


## Speed and memory

| corpus | reads | clusters | tool | peak RSS | time |
|---|---|---|---|---|---|
| SIRV real, full | 1 300 066 | 579 | python | 3.55 GB | 491 s |
| | | | isONclust3 | 2.49 GB | 243 s |
| | | | **port `--t 1`** | **0.74 GB** | **73 s** |
| | | | **port `--t 8`** | **1.34 GB** | **40 s** |
| Drosophila ONT | 1 000 000 | 84 318 
| | | | isONclust3 | 2.98 GB | 233 s |
| | | | **port `--t 1`** | **1.40 GB** | **726 s** |
| | | | **port `--t 8`** | **2.95 GB** | **205 s** |
| SIRV PacBio | 17 633 | 151 | isONclust3 | 0.12 GB | 1.7 s |
| | | | **port `--t 1`** | **0.75 GB** | **9.8 s** |
| | | | **port `--t 8`** | **0.60 GB** | **3.9 s** |

The Python implementation was not run on Drosophila at 1M reads (takes too long).

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
corpora. The `--t > 1` is slightly different from `--t 1`, as mentioend in the paper.

**Downstream isoform reconstruction on the SIRV PacBio clusters.** isONform run
over each tool's output (`bench/downstream.sh`), isoforms scored against the SIRV
reference:

| | isoforms | matching | recall | precision | F1 | isONform runtime |
|---|---|---|---|---|---|---|
| isONclust1 | 70 | 35 | 47.1% | 50.0% | 0.485 | 386 s |
| isONclust3 | 61 | 27 | 38.2% | 44.3% | 0.410 | 2053 s |

isONform's cost grows with cluster size, and the runtimes differ by ~1670 s.


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
