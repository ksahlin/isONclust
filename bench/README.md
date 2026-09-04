# bench — the equivalence harness

Byte-identity with the Python reference is the port's acceptance criterion. This directory is what
checks it. The reasoning behind all of it is in `../PORTING.md`; this is the operating manual.

## Setup, once

```bash
bench/setup_reference_env.sh                       # python 3.12 -> conda env isonclust-ref
bench/setup_reference_env.sh 3.11 isonclust-ref-311 # and 3.11, to reproduce Finding 1
```

parasail has no macOS arm64 wheel and is on no conda channel, so this builds it from source and works
around two undocumented traps in its `setup.py`. Read the comment at the top of the script before
changing it.

The reference environment is **parasail and pysam, and nothing else**. `spoa` used to be needed by
`--consensus`, which is now out of scope (`PORTING.md`, *Scope*), so the port needs no POA at all and
nothing here looks for it.

## Running

```bash
bench/equivalence.sh env      # is the reference usable? is sum() compensated?
bench/equivalence.sh seeds    # does the reference agree with itself? RUN THIS FIRST
bench/equivalence.sh cli record
bench/equivalence.sh record
bench/equivalence.sh verify   # needs rust/target/release/isONclust
bench/equivalence.sh dropped  # the port must refuse the out-of-scope flags
```

Override with `REF_PYTHON`, `PORT_BIN`, `CORPUS`, `GOLDEN`.

## Layout

| Path | What |
| --- | --- |
| `equivalence.sh` | the harness |
| `cases.tsv` | the case matrix. **TAB-separated — do not let an editor expand tabs** |
| `diffsummary.py` | column-aware diff; a line diff is useless on files carrying whole reads |
| `setup_reference_env.sh` | builds the pinned reference environment |
| `env/resolved-*.txt` | exact versions, per platform and interpreter |
| `golden/manifest.tsv` | per-file sha256 for every case, plus the environment fingerprint |
| `golden/cli/` | CLI contract: stdout, stderr and exit code per case, verbatim |
| `golden/sample/` | the default case's small outputs, kept readable by eye |

## Two things worth knowing before you trust a green run

**The goldens are hashes, not files.** Verbatim they are 318 MB. On a mismatch `verify` re-runs the
reference for that one case (~2 s) and produces a real diff, which is strictly better than storing
the files.

**A harness that has never failed has not been tested.** Both of this one's own bugs reported passes
(`PORTING.md`, *Finding 10*). Three gates are demonstrated to fail when they should, and should be
re-demonstrated after any change here:

```bash
# 1. `seeds` must FAIL on Python 3.11 -- the reference genuinely is seed-dependent there
REF_PYTHON=~/miniforge3/envs/isonclust-ref-311/bin/python bench/equivalence.sh seeds

# 2. `verify` must catch a single changed digit in one field of 1240 lines.
#    Point PORT_BIN at a stand-in that runs the reference and perturbs one error_rate.

# 3. `dropped` must distinguish three cases. Verified with three stand-ins:
#      refuses and names the flag  -> pass
#      accepts and ignores         -> fail ("accepted (exit 0)")
#      refuses without naming it   -> fail ("does not name it")
```

The third matters because "unrecognised argument" is the easy thing for a port to emit and is not
actionable for someone whose pipeline passes `--consensus`.

## Corpora

`bench/corpora.tsv` is the registry. `CORPUS` takes either a path or a name from it:

```bash
CORPUS=smoke          bench/equivalence.sh record   # the committed fixture
CORPUS=sirv_real_10k  bench/equivalence.sh record   # real ONT, 10k reads
CORPUS=droso_20k      bench/equivalence.sh record   # the one that reaches Finding 4
```

Only `smoke` is committed — `../test/sirv_sim_120.fastq`, 356 KB, 120 simulated SIRV reads at 7%
error spanning 54 of 68 transcripts, with the source transcript in every header. Everything else is
real data under `$ISONCLUST_DATA` (default `~/data/lrRNA-seq`), outside the repository.

**Goldens are corpus-specific.** `manifest.tsv` records the corpus sha256 in its header; changing
`CORPUS` means re-recording. Keep separate `GOLDEN` directories if you want more than one at a time.

### The fixture was replaced, and why that matters

`sample_alz_2k.fastq` (2500 PacBio CCS reads) was the corpus until it was measured. At `--k 15` it
produced **byte-identical output for `--w 15` and `--w 50`** — across every output file — despite a
19× difference in minimizer density and all 225 probability-table entries differing. A port with a
broken minimizer window would have passed all 27 cases.

| corpus | `--w 15` vs `--w 50` at `--k 15` | distinct results across the 24 `main` cases |
| --- | --- | --- |
| `sample_alz_2k` (old, 2500 reads) | byte-identical | 8 of 24 |
| `smoke` (new, 120 reads) | 70 vs 76 clusters | 15 of 24 |
| `sirv_real_10k` | 35 vs 763 clusters | not yet swept |

The smoke fixture is 20× smaller than what it replaced *and* discriminates better. It is still a
smoke test: `--q 0` and `--q 15` collapse onto the default on it, because simulated reads at a fixed
error rate have no quality spread. Real corpora are where the `--q` and threshold sweeps get their
power, and `droso_20k` is the only one that reaches the empty-string-minimizer path at all.

See `PORTING.md`, *Finding 4* and *Finding 5*.
