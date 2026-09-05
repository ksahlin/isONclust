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


def dump_mapping(args, out):
    """Record every get_best_cluster call the real driver makes.

    `cluster.get_best_cluster` is wrapped, not reimplemented: the sweep is
    stateful -- the minimizer database grows as reads become representatives --
    so the only faithful way to capture its inputs is to let the reference run
    and watch. This is PORTING.md method point 3, dumping from the live driver.

    Format, per call:

        CALL  <read_cl_id> <compressed_seq_len> <n_minimizers> <error_rate_read>
        CAND  <cl_id> <error_rate> <indices,...> <positions,...> <acc>
        ...one CAND per candidate, in the reference's dict order...
        RES   <best_cluster_id> <nr_shared_kmers> <mapped_ratio>

    Floats are written with repr() so the replay reads back the identical double.
    """
    import argparse as _argparse
    from modules import p_minimizers_shared

    real = cluster.get_best_cluster
    state = {"n": 0, "stop": False}

    def wrapper(read_cl_id, compressed_seq_len, hit_clusters_ids,
                hit_clusters_hit_positions, minimizers, nummber_of_minimizers,
                hit_clusters_hit_index, representatives, p_emp_probs, a):
        res = real(read_cl_id, compressed_seq_len, hit_clusters_ids,
                   hit_clusters_hit_positions, minimizers, nummber_of_minimizers,
                   hit_clusters_hit_index, representatives, p_emp_probs, a)
        if not state["stop"]:
            out.write("CALL\t{0}\t{1}\t{2}\t{3!r}\n".format(
                read_cl_id, compressed_seq_len, nummber_of_minimizers,
                representatives[read_cl_id][6]))
            for cl_id in hit_clusters_hit_positions:
                out.write("CAND\t{0}\t{1!r}\t{2}\t{3}\t{4}\n".format(
                    cl_id, representatives[cl_id][6],
                    ",".join(str(x) for x in hit_clusters_hit_index[cl_id]),
                    ",".join(str(x) for x in hit_clusters_hit_positions[cl_id]),
                    representatives[cl_id][2]))
            out.write("RES\t{0}\t{1}\t{2!r}\n".format(res[0], res[1], res[2]))
            state["n"] += 1
            if args.max_calls and state["n"] >= args.max_calls:
                state["stop"] = True
        return res

    cluster.get_best_cluster = wrapper
    try:
        p_min_shared = p_minimizers_shared.read_empirical_p()
        p_emp_probs = {}
        for k, w, p, e1, e2 in p_min_shared:
            if int(k) == args.k and abs(int(w) - args.w) <= 2:
                p_emp_probs[(float(e1), float(e2))] = float(p)
                p_emp_probs[(float(e2), float(e1))] = float(p)

        read_array = [(i, 0, acc, seq, qual, float(acc.split("_")[-1]))
                      for i, (acc, (seq, qual))
                      in enumerate(help_functions.readfq(open(args.sorted_fastq)))]
        clusters, representatives = {}, {}
        for i, b_i, acc, seq, qual, score in read_array:
            clusters[i] = [acc]
            representatives[i] = (i, b_i, acc, seq, qual, score)

        a = _argparse.Namespace(
            k=args.k, w=args.w, min_shared=args.min_shared,
            mapped_threshold=args.mapped_threshold,
            aligned_threshold=args.aligned_threshold,
            min_fraction=args.min_fraction,
            min_prob_no_hits=args.min_prob_no_hits,
            print_output=0 or 10 ** 9)
        cluster.reads_to_clusters(clusters, representatives, read_array,
                                  p_emp_probs, {}, 1, a)
    finally:
        cluster.get_best_cluster = real


def dump_parasail(args, out):
    """Record every parasail alignment the real driver performs.

    `cluster.parasail_block_alignment` is wrapped, so what is captured is what
    isONclust actually asks parasail for -- match 2, mismatch -2, gap_ext 1, and
    an opening penalty chosen per comparison from the two reads' summed error
    rates. isONform verified its parasail port against isONcorrect's parameters
    (match 4, mismatch -8, open 12); different penalties can reach different
    tie-breaking paths, so this checks isONclust's own call site.

    Format, one record per call:

        PARA\t<opening_penalty>\t<k>\t<match_id>\t<s1>\t<s2>\t<cigar>\t<alignment_ratio>
    """
    import argparse as _argparse
    from modules import p_minimizers_shared

    real = cluster.parasail_block_alignment
    state = {"n": 0}

    def wrapper(s1, s2, k, match_id, match_score=2, mismatch_penalty=-2,
                opening_penalty=5, gap_ext=1):
        # Recompute the cigar the same way the reference does, so the dump
        # carries it without changing what the driver receives.
        import parasail as _p
        m = _p.matrix_create("ACGT", match_score, mismatch_penalty)
        r = _p.sg_trace_scan_16(s1, s2, opening_penalty, gap_ext, m)
        if r.saturated:
            r = _p.sg_trace_scan_32(s1, s2, opening_penalty, gap_ext, m)
        cig = str(r.cigar.decode, "utf-8")
        res = real(s1, s2, k, match_id, match_score, mismatch_penalty,
                   opening_penalty, gap_ext)
        ratio = res[2][2]
        if not (args.max_calls and state["n"] >= args.max_calls):
            out.write("PARA\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6!r}\n".format(
                opening_penalty, k, match_id, s1, s2, cig, ratio))
            state["n"] += 1
        return res

    cluster.parasail_block_alignment = wrapper
    try:
        p_min_shared = p_minimizers_shared.read_empirical_p()
        p_emp_probs = {}
        for k, w, p, e1, e2 in p_min_shared:
            if int(k) == args.k and abs(int(w) - args.w) <= 2:
                p_emp_probs[(float(e1), float(e2))] = float(p)
                p_emp_probs[(float(e2), float(e1))] = float(p)
        read_array = [(i, 0, acc, seq, qual, float(acc.split("_")[-1]))
                      for i, (acc, (seq, qual))
                      in enumerate(help_functions.readfq(open(args.sorted_fastq)))]
        clusters, representatives = {}, {}
        for i, b_i, acc, seq, qual, score in read_array:
            clusters[i] = [acc]
            representatives[i] = (i, b_i, acc, seq, qual, score)
        a = _argparse.Namespace(
            k=args.k, w=args.w, min_shared=args.min_shared,
            mapped_threshold=args.mapped_threshold,
            aligned_threshold=args.aligned_threshold,
            min_fraction=args.min_fraction,
            min_prob_no_hits=args.min_prob_no_hits,
            print_output=10 ** 9)
        cluster.reads_to_clusters(clusters, representatives, read_array,
                                  p_emp_probs, {}, 1, a)
    finally:
        cluster.parasail_block_alignment = real


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", required=True, choices=["minimizers", "mapping", "parasail"])
    ap.add_argument("--sorted-fastq", required=True,
                    help="the reference's sorted.fastq -- NOT the raw input")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--w", type=int, required=True)
    ap.add_argument("--q", type=float, default=7.0)
    ap.add_argument("--min_shared", type=int, default=5)
    ap.add_argument("--mapped_threshold", type=float, default=0.7)
    ap.add_argument("--aligned_threshold", type=float, default=0.4)
    ap.add_argument("--min_fraction", type=float, default=0.8)
    ap.add_argument("--min_prob_no_hits", type=float, default=0.1)
    ap.add_argument("--max-calls", type=int, default=0,
                    help="stop after this many recorded calls (0 = no limit)")
    ap.add_argument("--out", default=None,
                    help="write the dump here instead of stdout. USE THIS: the "
                         "reference prints progress to stdout from inside "
                         "reads_to_clusters, so a stdout dump is interleaved "
                         "with 5 stray lines like 'Saved: 0 iterations.'")
    args = ap.parse_args()
    out = open(args.out, "w") if args.out else sys.stdout
    try:
        if args.stage == "minimizers":
            dump_minimizers(args.sorted_fastq, args.k, args.w, out)
        elif args.stage == "mapping":
            dump_mapping(args, out)
        elif args.stage == "parasail":
            dump_parasail(args, out)
    finally:
        if args.out:
            out.close()


if __name__ == "__main__":
    main()
