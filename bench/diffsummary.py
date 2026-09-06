#!/usr/bin/env python3
"""Summarise the difference between two of the tool's output files.

A whole-line diff is useless on these files: final_cluster_origins.tsv carries
the full read sequence and quality string in columns 3 and 4, so one differing
float in column 6 prints four kilobytes. This reports WHICH COLUMN moved, and by
how much, which is the part that tells you what broke.
"""
import sys

COLS = {
    "final_clusters.tsv": ["cluster_id", "read_acc"],
    "final_cluster_origins.tsv": ["cluster_id", "acc", "seq", "qual", "score", "error_rate"],
    "pre_clusters.csv": ["cluster_id", "read_acc"],
    "cluster_origins.csv": ["read_cl_id", "acc", "seq", "qual", "score", "error_rate"],
}


def trunc(s, n=28):
    return s if len(s) <= n else s[: n - 3] + "..."


def main(a_path, b_path, label=None, limit=6):
    name = label or a_path.rsplit("/", 1)[-1]
    names = COLS.get(name)
    a = open(a_path, encoding="utf-8", errors="replace").read().splitlines()
    b = open(b_path, encoding="utf-8", errors="replace").read().splitlines()

    if len(a) != len(b):
        print(f"        line count differs: {len(a)} vs {len(b)}")

    percol, shown, ndiff = {}, 0, 0
    for i, (la, lb) in enumerate(zip(a, b), 1):
        if la == lb:
            continue
        ndiff += 1
        fa, fb = la.split("\t"), lb.split("\t")
        if len(fa) != len(fb):
            percol.setdefault("<field count>", 0)
            percol["<field count>"] += 1
            continue
        for j, (x, y) in enumerate(zip(fa, fb)):
            if x == y:
                continue
            col = names[j] if names and j < len(names) else f"col{j + 1}"
            percol[col] = percol.get(col, 0) + 1
            if shown < limit:
                extra = ""
                # A float that differs only in its last digits is a summation
                # order artefact, not a logic error. Saying so saves an hour.
                try:
                    fx, fy = float(x), float(y)
                    if fx != fy:
                        rel = abs(fx - fy) / max(abs(fx), abs(fy), 1e-300)
                        extra = f"   rel={rel:.2e}"
                        if rel < 1e-12:
                            extra += "  (float rounding, not a logic difference)"
                except ValueError:
                    pass
                print(f"        line {i} [{col}]: {trunc(x)!r} vs {trunc(y)!r}{extra}")
                shown += 1

    print(f"        {ndiff} differing lines; columns touched: "
          + ", ".join(f"{k}x{v}" for k, v in sorted(percol.items(), key=lambda kv: -kv[1])))


if __name__ == "__main__":
    main(*sys.argv[1:])
