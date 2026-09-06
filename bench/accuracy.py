#!/usr/bin/env python3
"""Score a clustering against per-read ground truth.

The metrics are the ones isONclust's own paper uses (`scripts/compute_cluster_quality.py`):
homogeneity, completeness, V-measure and the adjusted Rand index. That script
takes truth from a sorted BAM, which means running an aligner; this takes it from
the read accession instead, which is exact and needs no alignment -- so it only
works on corpora whose reads carry their source, and that is the point of having
the simulated SIRV corpora.

They are implemented here rather than imported from scikit-learn so the harness
has no heavyweight dependency. `tests` validates them against scikit-learn when
it is installed, because hand-rolled information-theoretic metrics are exactly
the kind of thing that looks right and is off by a normalisation.

    bench/accuracy.py --clusters OUT/final_clusters.tsv --truth-from accession

Truth extraction:
  accession   `read_17_from_SIRV612` -> `SIRV612`   (simulated SIRV)
  suffix      everything after the last `_`
"""
import argparse
import math
import re
import sys
from collections import Counter, defaultdict


def entropy(counts, n):
    """H = -sum p log p, natural log, matching scikit-learn."""
    h = 0.0
    for c in counts:
        if c > 0:
            p = c / n
            h -= p * math.log(p)
    return h


def homogeneity_completeness_v(labels_true, labels_pred):
    """Rosenberg & Hirschberg (2007), as scikit-learn computes them."""
    n = len(labels_true)
    if n == 0:
        return 0.0, 0.0, 0.0
    ct = defaultdict(Counter)
    for t, p in zip(labels_true, labels_pred):
        ct[t][p] += 1
    true_counts = Counter(labels_true)
    pred_counts = Counter(labels_pred)

    h_c = entropy(true_counts.values(), n)
    h_k = entropy(pred_counts.values(), n)

    # H(C|K) = -sum_ck (n_ck/n) log(n_ck/n_k)
    h_ck = 0.0
    for t, row in ct.items():
        for p, n_ck in row.items():
            h_ck -= (n_ck / n) * math.log(n_ck / pred_counts[p])
    # H(K|C)
    h_kc = 0.0
    for t, row in ct.items():
        for p, n_ck in row.items():
            h_kc -= (n_ck / n) * math.log(n_ck / true_counts[t])

    # scikit-learn's convention: a degenerate entropy gives a score of 1.0
    homogeneity = 1.0 if h_c == 0.0 else 1.0 - h_ck / h_c
    completeness = 1.0 if h_k == 0.0 else 1.0 - h_kc / h_k
    if homogeneity + completeness == 0.0:
        v = 0.0
    else:
        v = 2.0 * homogeneity * completeness / (homogeneity + completeness)
    return homogeneity, completeness, v


def adjusted_rand_index(labels_true, labels_pred):
    n = len(labels_true)
    if n < 2:
        return 1.0
    n_classes = len(set(labels_true))
    n_clusters = len(set(labels_pred))
    # scikit-learn's special cases: no split at all, or every point its own
    # cluster, on BOTH sides, are perfect matches and score 1.0 rather than
    # falling out of the formula as 0/0. Found by the validation below, not by
    # reading the source first.
    if n_classes == n_clusters == 1 or n_classes == n_clusters == n:
        return 1.0
    ct = defaultdict(Counter)
    for t, p in zip(labels_true, labels_pred):
        ct[t][p] += 1
    comb2 = lambda x: x * (x - 1) / 2
    index = sum(comb2(v) for row in ct.values() for v in row.values())
    a = sum(comb2(v) for v in Counter(labels_true).values())
    b = sum(comb2(v) for v in Counter(labels_pred).values())
    total = comb2(n)
    expected = a * b / total
    maximum = (a + b) / 2
    if maximum == expected:
        return 0.0
    return (index - expected) / (maximum - expected)


TRUTH_PATTERNS = {
    # read_17_from_SIRV612 -> SIRV612
    "accession": re.compile(r"_from_(\S+?)(?:_|$)"),
}


def truth_of(acc, how):
    if how == "accession":
        m = TRUTH_PATTERNS["accession"].search(acc)
        if m:
            return m.group(1)
        return None
    if how == "suffix":
        return acc.rsplit("_", 1)[-1]
    raise SystemExit(f"unknown --truth-from {how!r}")


def read_clusters(path):
    """cluster_id -> [read accessions]. Both tools emit `<id>\\t<acc>`."""
    out = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            out.append((f[0], f[1]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--truth-file", default=None,
                    help="`<accession>\\t<class>` per line, from bench/make_truth.sh. "
                         "REQUIRED for real corpora, and the only honest option for "
                         "cross-tool comparison -- see that script's header on why the "
                         "simulated corpora cannot be used for it.")
    ap.add_argument("--truth-from", default="accession", choices=["accession", "suffix"])
    ap.add_argument("--label", default="")
    ap.add_argument("--tsv", action="store_true", help="one tab-separated line, for tables")
    args = ap.parse_args()

    pairs = read_clusters(args.clusters)

    truth = None
    if args.truth_file:
        truth = {}
        with open(args.truth_file) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2:
                    truth[f[0]] = f[1]

    labels_true, labels_pred, unknown = [], [], 0
    for cl, acc in pairs:
        t = truth.get(acc) if truth is not None else truth_of(acc, args.truth_from)
        if t is None:
            unknown += 1
            continue
        labels_true.append(t)
        labels_pred.append(cl)

    if not labels_true:
        where = args.truth_file or f"--truth-from {args.truth_from}"
        msg = (f"no read in {args.clusters} matched any truth from {where}; "
               f"{unknown} accessions did not match")
        if args.tsv:
            # A corpus the truth does not cover is normal when one benchmark run
            # spans several corpora. Report the cluster shape, which needs no
            # truth, and dashes for the metrics that do.
            sizes = Counter(cl for cl, _ in pairs)
            top10 = [v for _, v in sizes.most_common(10)]
            print("\t".join([args.label, "0", "0", str(len(sizes)),
                             str(sum(1 for v in sizes.values() if v == 1)),
                             str(max(sizes.values()) if sizes else 0),
                             "-", "-", "-", "-",
                             ",".join(str(x) for x in top10)]))
            print(f"# {msg}", file=sys.stderr)
            return
        raise SystemExit(msg + ". Check that the clustering and the truth use the "
                         "same accession form.")

    h, c, v = homogeneity_completeness_v(labels_true, labels_pred)
    ari = adjusted_rand_index(labels_true, labels_pred)

    # Cluster sizes are reported over EVERY cluster, not only the scored reads,
    # because a cluster made entirely of reads without truth is still a cluster
    # the tool produced and still costs the downstream stage.
    all_sizes = Counter(cl for cl, _ in pairs)
    n_clusters = len(all_sizes)
    singletons = sum(1 for v in all_sizes.values() if v == 1)
    largest = max(all_sizes.values())
    top10 = [v for _, v in all_sizes.most_common(10)]
    n_true = len(set(labels_true))

    if args.tsv:
        print("\t".join(str(x) for x in [
            args.label, len(labels_true), n_true, n_clusters, singletons, largest,
            f"{h:.6f}", f"{c:.6f}", f"{v:.6f}", f"{ari:.6f}",
            ",".join(str(x) for x in top10)]))
    else:
        print(f"  {args.label}")
        print(f"    reads scored      {len(labels_true)}  (unrecoverable truth: {unknown})")
        print(f"    true classes      {n_true}")
        print(f"    clusters          {n_clusters}  ({singletons} singletons, largest {largest})")
        print(f"    top 10 sizes      {', '.join(str(x) for x in top10)}")
        print(f"    homogeneity       {h:.4f}")
        print(f"    completeness      {c:.4f}")
        print(f"    V-measure         {v:.4f}")
        print(f"    adjusted Rand     {ari:.4f}")


if __name__ == "__main__":
    main()
