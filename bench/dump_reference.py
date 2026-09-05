#!/usr/bin/env python3
"""Dump a stage's inputs and outputs from the reference, without modifying it.

End-to-end equivalence says *that* the port is wrong, never *where*. Stages
whose output reaches a file can be diffed directly (`equivalence.sh stage sort`);
stages whose output stays in memory cannot, and this is how those get checked.

`get_kmer_minimizers` is the first such stage, and PORTING.md's Finding 5 is the
argument that it needs one: on the retired corpus, `--w 15` and `--w 50` produced
byte-identical output despite a 19x difference in minimizer density, so a port
with a broken window would have passed every end-to-end case.

The reference is imported, not copied or edited. If it changes, this changes
with it.

    bench/dump_reference.py --stage minimizers --sorted-fastq OUT/sorted.fastq \
        --k 13 --w 20 > minimizers.tsv

Output format, one line per minimizer, in the order the reference produces them:

    <read index>\t<position>\t<minimizer>

Read index is the 0-based position in sorted.fastq, which is the order
`reads_to_clusters` iterates. Reads the reference skips (homopolymer-compressed
length < k) emit a single line with position -1 and an empty minimizer, so the
skip itself is part of the contract rather than an absence.
"""
import argparse
import itertools
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from modules import help_functions  # noqa: E402
from modules import cluster  # noqa: E402


def dump_minimizers(sorted_fastq, k, w, out):
    for idx, (acc, (seq, qual)) in enumerate(help_functions.readfq(open(sorted_fastq))):
        # Exactly what reads_to_clusters does before calling the function.
        seq_hpol_comp = "".join(ch for ch, _ in itertools.groupby(seq))
        if len(seq_hpol_comp) < k:
            out.write("{0}\t-1\t\n".format(idx))
            continue
        for m, pos in cluster.get_kmer_minimizers(seq_hpol_comp, k, w):
            out.write("{0}\t{1}\t{2}\n".format(idx, pos, m))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", required=True, choices=["minimizers"])
    ap.add_argument("--sorted-fastq", required=True,
                    help="the reference's sorted.fastq -- NOT the raw input")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--w", type=int, required=True)
    args = ap.parse_args()
    if args.stage == "minimizers":
        dump_minimizers(args.sorted_fastq, args.k, args.w, sys.stdout)


if __name__ == "__main__":
    main()
