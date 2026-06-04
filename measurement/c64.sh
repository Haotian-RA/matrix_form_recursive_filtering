#!/bin/bash
# Compile and run PMCTest in 64 bit mode using C++
# (c) 2012 by Agner Fog. GNU General Public License www.gnu.org/licenses

# Auto-detect compiler
if command -v clang++ &>/dev/null; then
    CXX="clang++"
elif command -v g++ &>/dev/null; then
    CXX="g++"
else
    echo "Error: no C++ compiler found"
    exit 1
fi

echo "Using compiler: $CXX"

# Compile A file if modified
if [ PMCTestA.cpp -nt a64.o ] ; then
$CXX -O2 -c -m64 -oa64.o PMCTestA.cpp
fi

# Compile B file and link
$CXX -std=c++20 -mavx2 -mfma -O2 a64.o -march=native -ffast-math -lpthread PMCTestB.cpp
if [ $? -ne 0 ] ; then exit ; fi
taskset -c 0 ./a.out
