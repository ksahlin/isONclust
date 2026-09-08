# isONclust - Clustering of long-read transcriptome reads into gene families

### isONclust has been re-implemented in Rust (2026-09-08) and produces identical output 2-7x faster on a fraction of the memory (see below).

isONclust clusters PacBio Iso-Seq or Oxford Nanopore reads so that each cluster
holds the reads from one gene. Output is a tsv file assigning each read to a
cluster ID. Detailed information is in the [paper](https://link.springer.com/chapter/10.1007/978-3-030-17083-7_14).

## Installation <a name="installation"></a>

It needs a Rust toolchain ([rustup.rs](https://rustup.rs)), `cmake` and `libclang`.

```
git clone https://github.com/ksahlin/isONclust.git
cd isONclust/rust
cargo build --release
```

That produces `target/release/isONclust`. Put it on your `PATH` and run it as
shown under [Running isONclust](#Running). Without `cmake` or `libclang`, add
`--no-default-features`: same output, about 3x slower.

The original python implementation is still available and is the reference this
port is checked against; see [INSTALL-python.md](INSTALL-python.md).

### Rust-port versions

The port is a drop-in replacement: same command line, same flags, same output
files, **byte for byte**.

One divergence: sequences are 2-bit packed, so a non-ACGT base will be converted
to `A` (a warning is emitted if it sees any).

Peak RSS and runtime at `--t 1`:

| corpus | reads | python | Rust port |
|---|---|---|---|
| SIRV ONT | 10 000 | 235 MB / 4.25 s | **39 MB / 0.84 s** |
| Drosophila ONT | 20 000 | 704 MB / 16.1 s | **328 MB / 5.83 s** |
| SIRV PacBio | 17 633 | 1092 MB / 23.6 s | **746 MB / 9.80 s** |
| SIRV ONT, full | 1 300 066 | 3.55 GB / 491 s | **0.74 GB / 72.9 s** |

Full comparison --- accuracy against gene-level truth, cluster-size
distributions, and a comparison with
[isONclust3](https://github.com/aljpetri/isONclust3) --- is in
[Port-benchmark.md](Port-benchmark.md).

### Running a test <a name="runtest"></a>

`test/sirv_sim_120.fastq` is a small dataset for checking an installation.

```
isONclust --ont --fastq test/sirv_sim_120.fastq --outfolder /tmp/isonclust_test
```

This finishes in under a second and writes `final_clusters.tsv`. 

## Input data <a name="Input_data"></a>

A fastq with ONT or PacBio reads. Reads should have barcodes removed --- with
[LIMA](https://lima.how/) for PacBio, or
[Pychopper](https://github.com/epi2me-labs/pychopper) for ONT.


## Running isONclust <a name="Running"></a>


### ONT reads

```
isONclust --ont --fastq [reads.fastq] --outfolder [/path/to/output] 
```

The argument `--ont`  means `--k 13 --w 20`. These arguments can be set manually without the `--ont` flag. Specify number of cores with `--t`. 


### Iso-Seq reads


```
isONclust --isoseq --fastq [reads.fastq] --outfolder [/path/to/output]
```

The argument `--isoseq` simply means `--k 15 --w 50`. These arguments can be set manually without the `--isoseq` flag. Specify number of cores with `--t`. 


## Outputs <a name="Outputs"></a>

#### Clustering information
The output consists of a tsv file `final_clusters.tsv` present in the specified output folder. In this file, the first column is the cluster ID and the second column is the read accession. For example:
```
0 read_X_acc
0 read_Y_acc
...
n read_Z_acc
```
if there are n reads there will be n rows. Some reads might be singletons. The rows are ordered with respect to the size of the cluster (largest first).

#### Cluster fastq files

You can obtain separate cluster fastq files from the clustering by running

```
isONclust write_fastq --clusters [/path/to/output/]final_clusters.tsv --fastq [reads.fastq] --outfolder [/path/to/fastq_output] --N 1
```


## Credits <a name="credits"></a>

Please cite [1] when using isONclust.

1. Kristoffer Sahlin, Paul Medvedev. De Novo Clustering of Long-Read Transcriptome Data Using a Greedy, Quality-Value Based Algorithm, Journal of Computational Biology 2020, 27:4, 472-484. [Link](https://www.liebertpub.com/doi/abs/10.1089/cmb.2019.0299).

Here is an open access version of the paper: [bioRxiv link](https://www.biorxiv.org/content/10.1101/463463v1).

#### Bib record

@article{sahlin2020a,
author = {Sahlin, Kristoffer and Medvedev, Paul},
title = {De Novo Clustering of Long-Read Transcriptome Data Using a Greedy, Quality Value-Based Algorithm},
journal = {Journal of Computational Biology},
volume = {27},
number = {4},
pages = {472-484},
year = {2020},
doi = {10.1089/cmb.2019.0299},
note ={PMID: 32181688},
URL = {https://doi.org/10.1089/cmb.2019.0299},
eprint = {https://doi.org/10.1089/cmb.2019.0299},
abstract = { Long-read sequencing of transcripts with Pacific Biosciences (PacBio) Iso-Seq and Oxford Nanopore Technologies has proven to be central to the study of complex isoform landscapes in many organisms. However, current de novo transcript reconstruction algorithms from long-read data are limited, leaving the potential of these technologies unfulfilled. A common bottleneck is the dearth of scalable and accurate algorithms for clustering long reads according to their gene family of origin. To address this challenge, we develop isONclust, a clustering algorithm that is greedy (to scale) and makes use of quality values (to handle variable error rates). We test isONclust on three simulated and five biological data sets, across a breadth of organisms, technologies, and read depths. Our results demonstrate that isONclust is a substantial improvement over previous approaches, both in terms of overall accuracy and/or scalability to large data sets. }
}

## Licence

GPL v3.0, see [LICENSE.txt](https://github.com/ksahlin/isONclust/blob/master/LICENCE.txt).


