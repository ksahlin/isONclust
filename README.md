isONclust
========

isONclust is a tool for clustering either PacBio Iso-Seq reads, or Oxford Nanopore reads into clusters, where each cluster represents all reads that came from a gene. Output is a tsv file with each read assigned to a cluster-ID. Detailed information is available in [preprint](https://www.biorxiv.org/content/early/2018/11/06/463463).  


isONclust is distributed as a python package supported on Linux / OSX with python v>=3.4 as of version 0.0.2 and above (due to updates in python's multiprocessing library). [![Build Status](https://travis-ci.org/ksahlin/isONclust.svg?branch=master)](https://travis-ci.org/ksahlin/isONclust).

Table of Contents
=================

  * [INSTALLATION](#INSTALLATION)
    * [Using conda](#Using-conda)
    * [Using pip](#Using-pip)
    * [Downloading source from GitHub](#Downloading-source-from-github)
    * [Dependencies](#Dependencies)
    * [Testing installation](#testing-installation)
  * [USAGE](#USAGE)
    * [Iso-Seq](#Iso-Seq)
    * [Oxford Nanopore](#Oxford-Nanopore)
    * [Output](#Output)
    * [Parameters](#Parameters)
  * [CREDITS](#CREDITS)
  * [LICENCE](#LICENCE)



INSTALLATION
----------------

### Using conda
Conda is the preferred way to install isONclust.

1. Create and activate a new environment called isonclust

```
conda create -n isonclust python=3 pip 
source activate isonclust
```

2. Install isONclust 

```
pip install isONclust
```
3. You should now have 'isONclust' installed; try it:
```
isONclust --help
```

Upon start/login to your server/computer you need to activate the conda environment "isonclust" to run isONclust as:
```
source activate isonclust
```

### Using pip 

To install isONclust, run:
```
pip install  isONclust
```
`pip` will install the dependencies automatically for you. `pip` is pythons official package installer and is included in most python versions. If you do not have `pip`, it can be easily installed [from here](https://pip.pypa.io/en/stable/installing/) and upgraded with `pip install --upgrade pip`. 


### Downloading source from GitHub

#### Dependencies

Make sure the below listed dependencies are installed (installation links below). Versions in parenthesis are suggested as IsoCon has not been tested with earlier versions of these libraries. However, IsoCon may also work with earliear versions of these libaries.
* [parasail](https://github.com/jeffdaily/parasail-python)
* [pysam](http://pysam.readthedocs.io/en/latest/installation.html) (>= v0.11)


With these dependencies installed. Run

```sh
git clone https://github.com/ksahlin/isONclust.git
cd isONclust
./isONclust
```

### Testing installation

You can verify successul installation by running isONclust on this [small dataset](https://github.com/ksahlin/isONclust/tree/master/test/sample_alz_2k.fastq). Simply download the test dataset and run:

```
isONclust --fastq [test/sample_alz_2k.fastq] --outfolder [output path]
```


USAGE
-------

IsONclust can be used with either Iso-Seq or ONT reads. It takes either a fastq file or ccs.bam file. 
 


### Iso-Seq

IsONclust works with full-lengh non-chimeric (_flnc_) reads that has quality values assigned to bases. The flnc reads with quality values can be generated as follows:

1. Make sure quality values is output when running the circular consensus calling step (CCS), by running `ccs` with the parameter `--polish`.
2. Run PacBio's Iso-Seq pipeline step 2 and 3 (primer removal and extraction of flnc reads) [isoseq3](https://github.com/PacificBiosciences/IsoSeq3/blob/master/README_v3.1.md).  

Flnc reads can be submitted as either a fastq file or bam file. A fastq file is created from a BAM by running _e.g_ `bamtools convert -format fastq -in flnc.bam -out flnc.fastq`. isONclust is called as follows

```
isONclust pipeline --isoseq --fastq <reads.fastq> --outfolder </path/to/output> 
```

isONclust also supports older versions of the isoseq3 pipeline by taking the `ccs.bam` file together with the `flnc.bam`. In this case, isONclust can be run as follows. 

<!--- If not, flnc reads can be generated as follows. Raw pacbio subreads needs to be proccesed with `ccs` with the command `--polish` (to get quality values), followed by `lima`, and `isoseq3 cluster` to get the flnc reads. The flnc file is generated at the very beginning of the `isoseq3 cluster` algorithm and it can be used once its created (no need to wait for isoseq3 to finish). See full documentation on generating flnc reads at [isoseq3](https://github.com/PacificBiosciences/IsoSeq3). After these three comands are run isONclust can be run as follows -->
```
isONclust --isoseq --ccs <ccs.bam> --flnc <flnc.bam> --outfolder </path/to/output> 
```
Where `<ccs.bam>` is the file generated from `ccs` and `<flnc.bam>` is the file generated from `isoseq3 cluster`. The argument `--isoseq` simply means `--k 15 --w 50`. These arguments can be set manually without the `--isoseq` flag. Specify number of cores with `--t`. 


### Oxford Nanopore
isONclust needs a fastq file generated by an Oxford Nanopore basecaller.

```
IsoCon pipeline --ont --fastq <reads.fastq> --outfolder </path/to/output> 
```
The argument `--ont` simply means `--k 13 --w 20`. These arguments can be set manually without the `--ont` flag. Specify number of cores with `--t`. 

## Output

### TSV

The output consists of a tsv file `final_clusters.tsv` present in the specified output folder. In this file, the first column is the cluster ID and the second column is the read accession. For example:
```
0 read_X_acc
0 read_Y_acc
...
n read_Z_acc
```
if there are n reads there will be n rows. Some reads might be singletons. The rows are ordered with respect to the size of the cluster (largest first). 

### Fastq

isONclust can also print separate fastq files for each cluster with more than N reads (N is a parameter to the program). After clustering, simply run
```
isONclust write_fastq --N [int] --fastq <reads.fastq> --clusters <path/to/final_clusters.tsv> --outfolder </path/to/output>
```
This will print out separate fastq files in `</path/to/output>` for all clusters with more than `[int]` reads. The names of the files are the cluster IDs assigned by isONclust, and matches the ID's found in the `final_clusters.tsv` file.


### Parameters

```
optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --fastq FASTQ         Path to consensus fastq file(s) (default: False)
  --flnc FLNC           The flnc reads generated by the isoseq3 algorithm (BAM
                        file) (default: False)
  --ccs CCS             Path to consensus BAM file(s) (default: False)
  --t NR_CORES          Number of cores allocated for clustering (default: 8)
  --ont                 Clustering of ONT transcript reads. (default: False)
  --isoseq              Clustering of PacBio Iso-Seq reads. (default: False)
  --k K                 Kmer size (default: 15)
  --w W                 Window size (default: 50)
  --min_shared MIN_SHARED
                        Minmum number of minimizers shared between read and
                        cluster (default: 5)
  --mapped_threshold MAPPED_THRESHOLD
                        Minmum mapped fraction of read to be included in
                        cluster. The density of minimizers to classify a
                        region as mapped depends on quality of the read.
                        (default: 0.7)
  --aligned_threshold ALIGNED_THRESHOLD
                        Minmum aligned fraction of read to be included in
                        cluster. Aligned identity depends on the quality of
                        the read. (default: 0.4)
  --min_fraction MIN_FRACTION
                        Minmum fraction of minimizers shared compared to best
                        hit, in order to continue mapping. (default: 0.8)
  --min_prob_no_hits MIN_PROB_NO_HITS
                        Minimum probability for i consecutive minimizers to be
                        different between read and representative and still
                        considered as mapped region, under assumption that
                        they come from the same transcript (depends on read
                        quality). (default: 0.1)
  --outfolder OUTFOLDER
                        A fasta file with transcripts that are shared between
                        samples and have perfect illumina support. (default:
                        None)
```

CREDITS
----------------

Please cite [1] when using IsoCon.

1. Kristoffer Sahlin, Paul Medvedev (2018) "De novo clustering of long-read transcriptome data using a greedy, quality-value based algorithm", bioRxiv [Link](https://www.biorxiv.org/content/early/2018/11/06/463463).

LICENCE
----------------

GPL v3.0, see [LICENSE.txt](https://github.com/ksahlin/IsoCon/blob/master/LICENCE.txt).

