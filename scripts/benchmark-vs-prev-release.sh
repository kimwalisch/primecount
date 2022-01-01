#!/bin/bash

# Benchmark the latest code vs. the previous release.
# Usage:
#   ./benchmark-vs-prev-release.sh
# Description:
#   Benchmarks the latest code vs the previous release.
#   Computes PrimePi(1e17) 10 times with each binary and sum
#   the elapsed seconds for each binary. The new code must
#   not be more than 3% slower.

# Exit if any error occurs
set -e

# Execute in base directory
if [[ "$(basename $(pwd))" = "scripts" ]]
then
    cd ..
fi

command -v bc >/dev/null 2>/dev/null
if [[ $? -ne 0 ]]
then
    echo "Error: GNU bc is not installed."
    exit 1
fi

rm -rf build-curr-release build-prev-release 

# Build the latest code
cmake -S . -B build-curr-release  -G "Unix Makefiles"
cmake --build build-curr-release -- -j4

# Checkout the previous release tag
git checkout $(git describe --tags --abbrev=0)
cmake -S . -B build-prev-release  -G "Unix Makefiles"
cmake --build build-prev-release -- -j4

# New code must not be more than 
# 3% slower than old code.
factor=1.03

# Test failure must be observed 3 times,
# we try to avoid false negatives.
for j in {1..3}
do
    total_seconds1=0
    total_seconds2=0

    for i in {1..10}
    do
        seconds1=$(build-prev-release/./primecount 1e17 --time | grep Seconds | cut -d' ' -f2)
        echo "Seconds old code: $seconds1"
        total_seconds1=$(echo "scale=3; $total_seconds1 + $seconds1" | bc -l)
        sleep 1
        seconds2=$(build-curr-release/./primecount 1e17 --time | grep Seconds | cut -d' ' -f2)
        echo "Seconds new code: $seconds2"
        echo ""
        total_seconds2=$(echo "scale=3; $total_seconds2 + $seconds2" | bc -l)
        sleep 1
    done

    echo "Total seconds old code: $total_seconds1"
    echo "Total seconds new code: $total_seconds2"
    echo ""

    limit=$(echo "scale=3; $total_seconds1 * $factor" | bc -l)
    new_code_is_fast=$(echo $total_seconds2'<='$limit | bc -l)

    if [ $new_code_is_fast -eq 1 ]
    then
        new_code_percent=$(echo "scale=1; 100 * $total_seconds1 / $total_seconds2" | bc -l)
        echo "Old code: 100.0%"
        echo "New code: $new_code_percent%"
        echo "New code successfully passed performance test!"
        exit 0
    fi
done

echo "New code is more than $factor times slower than previous release!"
exit 1
