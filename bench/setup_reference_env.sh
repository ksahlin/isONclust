#!/usr/bin/env bash
#
# Build a pinned conda environment that can run the isONclust Python reference.
#
#   bench/setup_reference_env.sh                 # default: python 3.12, env isonclust-ref
#   bench/setup_reference_env.sh 3.11 isonclust-ref-311
#
# Why this is not just `pip install -r requirements.txt`
# -----------------------------------------------------
# parasail publishes NO wheel for macOS arm64 (1.3.4 ships macosx_10_9_x86_64
# only) and it is not on conda-forge or bioconda under any name. On Apple
# Silicon it must be built from source, and its setup.py fights you twice:
#
#   1. It PREPENDS /usr/bin to PATH before probing for build tools, so macOS's
#      m4 1.4.6 always beats a newer m4 on your PATH. When it decides m4 is too
#      old it downloads m4-1.4.17 from ftp.gnu.org and builds it -- which fails.
#      Setting $M4 is the way out, because autoreconf honours it.
#   2. On Darwin it requires the tools named `glibtoolize` and `glibtool`
#      (the Homebrew spelling). conda-forge's libtool installs them as
#      `libtoolize` and `libtool`, so the probe fails on a perfectly good
#      libtool. Two symlinks fix it.
#
# Neither failure names its cause: the first reports "autoreconf -fi failed"
# after a detour through building m4, and the second reports nothing at all.
#
# The interpreter version is not a detail -- see PORTING.md, "sum() is
# compensated from 3.12 and naive before it". It changes whether the reference
# is deterministic. Build both if you want to reproduce that measurement.
set -euo pipefail

PYVER="${1:-3.12}"
ENVNAME="${2:-isonclust-ref}"
PARASAIL_VERSION="1.3.4"

CONDA="${CONDA_EXE:-$(command -v conda || true)}"
if [[ -z "$CONDA" ]]; then
  for c in "$HOME/miniforge3/bin/conda" "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda"; do
    [[ -x "$c" ]] && CONDA="$c" && break
  done
fi
[[ -n "$CONDA" ]] || { echo "error: conda not found. Install miniforge." >&2; exit 1; }

CONDA_ROOT="$(dirname "$(dirname "$CONDA")")"
ENVDIR="$CONDA_ROOT/envs/$ENVNAME"

echo "==> creating env '$ENVNAME' (python $PYVER) with the autotools parasail needs"
"$CONDA" create -y -n "$ENVNAME" -c conda-forge \
  "python=$PYVER" autoconf automake libtool pkg-config make m4 cxx-compiler \
  2>&1 | tail -3

echo "==> shimming glibtoolize/glibtool (parasail's Darwin probe uses the Homebrew names)"
SHIM="$ENVDIR/shim"
mkdir -p "$SHIM"
ln -sf "$ENVDIR/bin/libtoolize" "$SHIM/glibtoolize"
ln -sf "$ENVDIR/bin/libtool"    "$SHIM/glibtool"
ln -sf "$ENVDIR/bin/m4"         "$SHIM/m4"

export PATH="$SHIM:$ENVDIR/bin:$PATH"
export M4="$ENVDIR/bin/m4"
export LIBTOOLIZE="$SHIM/glibtoolize"
export LIBTOOL="$SHIM/glibtool"

echo "==> installing pysam (wheel) and parasail==$PARASAIL_VERSION (source build on arm64)"
"$ENVDIR/bin/pip" install -q --upgrade pip
"$ENVDIR/bin/pip" install -q pysam
"$ENVDIR/bin/pip" install "parasail==$PARASAIL_VERSION" --no-cache-dir --no-build-isolation 2>&1 | tail -3

echo "==> verifying"
"$ENVDIR/bin/python" - <<'PY'
import sys, parasail, pysam
m = parasail.matrix_create("ACGT", 2, -2)
r = parasail.sg_trace_scan_16("ACGTACGTAA", "ACGTTCGTAA", 5, 1, m)
assert r.score == 16 and str(r.cigar.decode, "utf-8") == "4=1X5=", "parasail behaves unexpectedly"
assert not r.saturated
print("    python  ", sys.version.split()[0])
print("    pysam   ", pysam.__version__)
print("    parasail  sg_trace_scan_16 OK (score 16, cigar 4=1X5=)")
# The interpreter's summation behaviour is part of the pinned contract.
xs = [0.1] * 10 + [1e17, -1e17]
print("    sum() is", "COMPENSATED (>=3.12)" if sum(xs) == sum(reversed(xs)) else "NAIVE (<=3.11)")
PY

echo
echo "==> done. Record the resolved versions:"
echo "    $ENVDIR/bin/pip freeze > bench/env/resolved-\$(uname -s)-\$(uname -m)-py$PYVER.txt"
echo "==> run the reference with:"
echo "    $ENVDIR/bin/python isONclust --ont --fastq test/sirv_sim_120.fastq --outfolder /tmp/out --t 1"
