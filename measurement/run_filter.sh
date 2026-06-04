#!/bin/bash
# run_filter.sh — Compile and run the IIR filter benchmark
#
# Usage:
#   ./run_filter.sh [options]
#
# Options:
#   --algo  scalar|bf|ph|cr        Algorithm reference (default: cr)
#   --arch  haswell|skylake|meteorlake  Architecture reference (default: haswell)
#   --order  2|4|8|16              Filter order (default: 16)
#   --blocksize  N                 Block count (default: 8)
#   --vector-size  SIZE            Input samples (default: 131072)
#   --iters  N                     Benchmark iterations (default: 10000)
#
# Examples:
#   ./run_filter.sh                                          # defaults: CR, order 16, N=8
#   ./run_filter.sh --algo ph --order 8 --blocksize 16       # PH, order 8, N=16
#   ./run_filter.sh --algo cr --order 16 --arch meteorlake   # CR on Meteor Lake

set -euo pipefail

# ─── Defaults ────────────────────────────────────────────────────────────────
ALGO="cr"
ARCH="haswell"
ORDER=16
BLOCKSIZE=8
VECSIZE=131072
ITERS=10000

# ─── Parse arguments ─────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case $1 in
        --algo)       ALGO="$2";      shift 2 ;;
        --arch)       ARCH="$2";      shift 2 ;;
        --order)      ORDER="$2";     shift 2 ;;
        --blocksize)  BLOCKSIZE="$2"; shift 2 ;;
        --vector-size) VECSIZE="$2";  shift 2 ;;
        --iters)      ITERS="$2";     shift 2 ;;
        --help|-h)
            head -20 "$0" | grep '^#' | sed 's/^# \?//'
            exit 0 ;;
        *)
            echo "Unknown option: $1 (try --help)"
            exit 1 ;;
    esac
done

# ─── Validate ────────────────────────────────────────────────────────────────
M=$((ORDER / 2))
if [[ "$M" -lt 1 || "$M" -gt 8 ]]; then
    echo "Error: --order must be 2, 4, 8, or 16"
    exit 1
fi

# ─── Enable perf_event if needed ─────────────────────────────────────────────
PARANOID=$(cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null || echo "unknown")
if [[ "$PARANOID" != "0" && "$PARANOID" != "-1" ]]; then
    echo "Setting perf_event_paranoid=0 (requires sudo)..."
    sudo sysctl -w kernel.perf_event_paranoid=0
    echo ""
fi

# ─── Auto-detect compiler ───────────────────────────────────────────────────
if command -v clang++ &>/dev/null; then
    CXX="clang++"
elif command -v g++ &>/dev/null; then
    CXX="g++"
else
    echo "Error: no C++ compiler found"
    exit 1
fi

# ─── Compile ─────────────────────────────────────────────────────────────────
SRC="$(cd "$(dirname "$0")" && pwd)/run_filter.cpp"
BIN="/tmp/run_filter_M${M}_N${BLOCKSIZE}"

echo "Compiling: M=$M (order $ORDER), N=$BLOCKSIZE, vector_size=$VECSIZE"
echo "  $CXX -> $BIN"

$CXX -std=c++20 -mavx2 -mfma -march=native -O2 -ffast-math -lpthread \
    -DM_VAL=$M -DN_VAL=$BLOCKSIZE -DVECTOR_SIZE=$VECSIZE -DITERS_VAL=$ITERS \
    "$SRC" -o "$BIN"

echo ""

# ─── Run ─────────────────────────────────────────────────────────────────────
taskset -c 0 "$BIN" --algo "$ALGO" --arch "$ARCH"
