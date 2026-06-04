#!/bin/bash
# Run all correctness tests (doctest-based).
# Usage: chmod +x test.sh && ./test.sh

set -uo pipefail

if command -v clang++ &>/dev/null; then
    CXX="clang++"
elif command -v g++ &>/dev/null; then
    CXX="g++"
else
    echo "Error: no C++ compiler found"
    exit 1
fi
CXXFLAGS="-std=c++20 -mavx2 -mfma -march=native -O2 -ffast-math"
BUILDDIR="/tmp/iir_tests"

mkdir -p "$BUILDDIR"

# All doctest test files (excludes benchmarks)
TESTS=(
    block_filtering.cpp
    fir_cores.cpp
    permute.cpp
    ph_decompos.cpp
    cyclic_reduction.cpp
    iir_cores.cpp
    filter.cpp
)

passed=0
failed=0
compile_failed=0
total=${#TESTS[@]}

echo "========================================"
echo "  Running $total correctness tests"
echo "========================================"
echo ""

for src in "${TESTS[@]}"; do
    name="${src%.cpp}"
    bin="$BUILDDIR/$name"

    printf "%-24s" "$name"

    # Compile
    if ! $CXX $CXXFLAGS "$src" -o "$bin" 2>/tmp/iir_compile_err.txt; then
        echo "COMPILE FAILED"
        cat /tmp/iir_compile_err.txt
        echo ""
        ((compile_failed++))
        continue
    fi

    # Run
    if output=$("$bin" 2>&1); then
        echo "PASSED"
        ((passed++))
    else
        echo "FAILED"
        echo "$output" | tail -20
        echo ""
        ((failed++))
    fi
done

echo ""
echo "========================================"
echo "  Results: $passed passed, $failed failed, $compile_failed compile errors (of $total)"
echo "========================================"

exit $((failed + compile_failed))
