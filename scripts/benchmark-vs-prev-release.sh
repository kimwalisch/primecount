#!/bin/bash

# Benchmark the latest code vs. the previous release.
# Usage:
#   ./benchmark-vs-prev-release.sh
# Description:
#   Benchmarks the latest code vs the previous release.
#   Computes PrimePi(1e17) 10 times with each binary and sum
#   the elapsed seconds for each binary. The new code must
#   not be more than 2% slower.

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
current_git_branch=$(git rev-parse --abbrev-ref HEAD)

# Build the latest code
cmake -S . -B build-curr-release  -G "Unix Makefiles"
cmake --build build-curr-release -- -j4

# Checkout the previous release tag
git checkout $(git describe --tags --abbrev=0)
cmake -S . -B build-prev-release  -G "Unix Makefiles"
cmake --build build-prev-release -- -j4

# Go back to initial branch
git checkout $current_git_branch

echo "=== Old code (previous release) ==="
build-prev-release/./primecount --version
echo ""

echo "=== New code ==="
build-curr-release/./primecount --version
echo ""

function benchmark_test_option {
    # New code must not be more than 
    # 5% slower than old code.
    factor=1.05

    # Test failure must be observed 3 times,
    # we try to avoid false negatives.
    for j in {1..3}
    do
        echo "=== Benchmark primecount --test ===="
        echo ""

        start=$(date +%s.%N)
        build-prev-release/./primecount -t4 --test >/dev/null
        end=$(date +%s.%N)
        seconds1=$(echo "scale=3; ($end - $start) / 1" | bc -l)
        echo "Seconds old code: $seconds1"
        sleep 1

        start=$(date +%s.%N)
        build-curr-release/./primecount -t4 --test >/dev/null
        end=$(date +%s.%N)
        seconds2=$(echo "scale=3; ($end - $start) / 1" | bc -l)
        echo "Seconds new code: $seconds2"
        sleep 1

        limit=$(echo "scale=3; $seconds1 * $factor" | bc -l)
        new_code_is_fast=$(echo $seconds2'<='$limit | bc -l)
        new_code_percent=$(echo "scale=1; 100 * $seconds2 / $seconds1" | bc -l)

        echo ""
        echo "Old code: 100.0%"
        echo "New code: $new_code_percent%"

        if [ $new_code_is_fast -eq 1 ]
        then
            echo "primecount --test performance test passed successfully!"
            return
        fi
    done

    echo "primecount --test is more than $factor times slower than previous release!"
    exit 1
}

function benchmark_pi_1e17 {
    # New code must not be more than 
    # 2% slower than old code.
    factor=1.02

    # Test failure must be observed 3 times,
    # we try to avoid false negatives.
    for j in {1..3}
    do
        echo "=== Benchmark PrimePi(1e17) ===="
        echo ""

        total_seconds1=0
        total_seconds2=0

        for i in {1..10}
        do
            seconds1=$(build-prev-release/./primecount 1e17 -t4 --time | grep Seconds | cut -d' ' -f2)
            echo "Seconds old code: $seconds1"
            total_seconds1=$(echo "scale=3; ($total_seconds1 + $seconds1) / 1" | bc -l)
            sleep 1
            seconds2=$(build-curr-release/./primecount 1e17 -t4 --time | grep Seconds | cut -d' ' -f2)
            echo "Seconds new code: $seconds2"
            echo ""
            total_seconds2=$(echo "scale=3; ($total_seconds2 + $seconds2) / 1" | bc -l)
            sleep 1
        done

        echo "Total seconds old code: $total_seconds1"
        echo "Total seconds new code: $total_seconds2"
        echo ""

        limit=$(echo "scale=3; $total_seconds1 * $factor" | bc -l)
        new_code_is_fast=$(echo $total_seconds2'<='$limit | bc -l)
        new_code_percent=$(echo "scale=1; 100 * $total_seconds2 / $total_seconds1" | bc -l)

        echo "Old code: 100.0%"
        echo "New code: $new_code_percent%"

        if [ $new_code_is_fast -eq 1 ]
        then
            echo "PrimePi(1e17) performance test passed successfully!"
            return
        fi
    done

    echo "PrimePi(1e17) is more than $factor times slower than previous release!"
    exit 1
}

# Execute benchmark functions
benchmark_test_option
echo ""
benchmark_pi_1e17
