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
expected to produce different clusters. Read the comparison as "which tool for
which data", not as a scoreboard — and see *Where isONclust3 struggles, and why
it is probably fixable* below before drawing conclusions from the SIRV rows.

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
a target and nothing should be tuned to improve them: they pick a different
winner on ONT than on PacBio, while the gene-level verdict is stable.

## Speed and memory

| corpus | preset | tool | `--t` | secs | peak MB | speedup |
|---|---|---|---|---|---|---|
| SIRV ONT, 10k | ont | python | 1 | 4.25 | 235 | — |
| | | **port** | 1 | **0.80** | **119** | **5.3x** |
| | | python | 8 | 1.33 | 206 | — |
| | | **port** | 8 | **0.25** | 177 | **5.3x** |
| | | isONclust3 | — | 0.65 | 32 | |
| SIRV PacBio, 17.6k | isoseq | python | 1 | 23.64 | 1092 | — |
| | | **port** | 1 | **9.53** | 1082 | **2.5x** |
| | | python | 8 | 7.06 | 454 | — |
| | | **port** | 8 | **3.68** | 1051 | **1.9x** |
| | | isONclust3 | — | 1.64 | 119 | |
| Drosophila ONT, 20k | ont | python | 1 | 16.08 | 704 | — |
| | | **port** | 1 | **6.29** | **525** | **2.6x** |
| | | python | 8 | 7.46 | 826 | — |
| | | **port** | 8 | **3.08** | **740** | **2.4x** |
| | | isONclust3 | — | 1.62 | 204 | |

**The port is 1.9–5.3x faster than the reference on every corpus and thread
count**, at equal or lower memory except on PacBio at `--t 8`, where eight
resident batches cost more than the reference's eight processes.

**Memory is where both isONclust versions lose to isONclust3**, by 2.6x on
Drosophila, 3.7x on SIRV ONT and 9.1x on SIRV PacBio. That is inherited from the
algorithm, not from the port: the reference holds every read's sequence *and*
quality string resident, and holds them more than once — the sorted file is read
back into a second array, and `reads_to_clusters` copies each entry again. The
port reproduces that structure because reproducing it exactly is the contract.
Packing nucleotides two bits each and caching one float per read instead of
retaining the quality string are the two obvious remedies; both are specified,
with their costs, under *Memory: the whole dataset is resident* in
[PORTING.md](PORTING.md), and neither has been implemented.

Most of that speed is one change: linking parasail's C library instead of using
the port's own exact scalar reimplementation of it. Alignment is 96–99.6% of
runtime, and the scalar version is 13–16x slower than the C library. Before that
change the port was *slower* than the reference on PacBio and Drosophila. Both
paths are exact; `--no-default-features` builds the pure-Rust one, which needs no
cmake or libclang and gives identical output more slowly.

## Accuracy — gene level

### SIRV ONT, 9998 reads with truth, 7 genes

| tool | `--t` | clusters | homogeneity | completeness | V | ARI | max cluster |
|---|---|---|---|---|---|---|---|
| python / **port** | 1 | 36 | 1.0000 | 0.5729 | 0.7285 | 0.5681 | 3835 |
| python / **port** | 8 | 30 | 1.0000 | 0.6472 | **0.7858** | **0.7338** | 4595 |
| isONclust3 | — | 66 | 1.0000 | 0.5005 | 0.6671 | 0.3119 | 1776 |

### SIRV PacBio, 14 783 reads with truth, 7 genes

| tool | `--t` | clusters | homogeneity | completeness | V | ARI | max cluster |
|---|---|---|---|---|---|---|---|
| python / **port** | 1 | 151 | 1.0000 | 0.6853 | 0.8132 | 0.7138 | 3327 |
| python / **port** | 8 | 110 | 1.0000 | 0.6933 | **0.8189** | **0.7259** | 3334 |
| isONclust3 | — | 163 | **0.6482** | 0.6503 | 0.6492 | 0.3578 | **8237** |

### Drosophila ONT, 20 000 reads, 3871 genes

| tool | `--t` | clusters | homogeneity | completeness | V | ARI | max cluster |
|---|---|---|---|---|---|---|---|
| python / **port** | 1 | 5679 | 0.9932 | 0.9826 | 0.9878 | 0.9285 | 459 |
| python / **port** | 8 | 5538 | 0.9914 | 0.9882 | **0.9898** | 0.9388 | 488 |
| isONclust3 | — | 6424 | 0.9937 | 0.9676 | 0.9805 | **0.9563** | 526 |

### Top 10 cluster sizes

| corpus | tool | sizes |
|---|---|---|
| SIRV ONT | isONclust1 `--t 8` | 4595, 944, 641, 567, 512, 423, 381, 329, 267, 209 |
| | isONclust3 | 1776, 1632, 1067, 655, 644, 489, 449, 423, 330, 329 |
| SIRV PacBio | isONclust1 `--t 8` | 3334, 3146, 1104, 1085, 1030, 866, 717, 573, 550, 355 |
| | isONclust3 | **8237**, 1347, 1021, 861, 805, 734, 571, 464, 440, 435 |
| Drosophila | isONclust1 `--t 8` | 488, 432, 313, 267, 188, 173, 160, 142, 140, 130 |
| | isONclust3 | 526, 328, 296, 277, 174, 171, 168, 147, 146, 134 |

## Reading it

**On Drosophila — transcriptome scale, 3871 genes — the two tools are close and
both are good.** V 0.990 against 0.981, and isONclust3 takes the ARI (0.956
against 0.939). This is the realistic case: thousands of genes at modest depth,
where the job is to avoid over-merging. Cluster shapes are similar (max 488
against 526). On this evidence there is no strong quality argument between them,
and isONclust3 is 2–5x faster again than the port.

**On SIRV, isONclust1 wins clearly, but SIRV is an unusual clustering problem.**
Seven genes and 10–18k reads means ~1500–2500 reads per gene, so the task is
almost entirely *completeness* — merge aggressively into very few, very deep
clusters. That is the axis isONclust3 does worst on. A 7-class truth is not
representative of a transcriptome, and the SIRV numbers should not be read as a
general verdict.

**`--t 8` is not just faster than `--t 1` but slightly more accurate**, on all
three corpora. The `--t > 1` path is a hierarchical batch-and-merge rather than a
parallelised single pass, and it appears to merge somewhat better. Unexplained,
and worth understanding rather than relying on.

## Where isONclust3 struggles, and why it is probably fixable

The one clearly bad result is **SIRV PacBio**, where isONclust3's homogeneity
falls to 0.6482 — its clusters mix genes — and its largest cluster is **8237
reads against isONclust1's 3327**. Everywhere else its homogeneity is ≥0.99.

That single cluster is expensive downstream. Running isONform over each tool's
PacBio clusters (`bench/downstream.sh`, isoforms scored against the SIRV
reference):

| | isoforms | matching | recall | precision | F1 | isONform runtime |
|---|---|---|---|---|---|---|
| **isONclust1** | 70 | 35 | **47.1%** | **50.0%** | **0.485** | **386 s** |
| isONclust3 | 61 | 27 | 38.2% | 44.3% | 0.410 | **2053 s** |

isONclust3 saves ~5 seconds of clustering and costs ~1670 seconds of isoform
reconstruction, because isONform's cost grows steeply with cluster size.

**This looks like a fixable gap rather than a fundamental one.** isONclust1's
homogeneity is a perfect 1.0000 on every SIRV row, and the mechanism is not
mysterious: when minimizer sharing is ambiguous it falls back to an *alignment*
before committing two reads to the same cluster, and that alignment is what stops
distinct genes being merged. isONclust3 has no such step. Adding an alignment
check — even a cheap one, gated to the cases where a merge would create a large
or low-identity cluster — would plausibly close most of this, and the parasail
binding described above shows the cost of exact affine alignment is far lower
than it is usually assumed to be. **This should be read as a suggestion for
isONclust3, not as a verdict against it.**

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
