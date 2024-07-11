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

function benchmark_test_option {
    # New code must not be more than 
    # 3% slower than old code.
    factor=1.03

    # Test failure must be observed 4 times,
    # we try to avoid false negatives.
    for j in {1..4}
    do
        echo ""
        echo "=== Benchmark primecount --test ===="
        echo ""

        start=$(date +%s.%N)
        build-prev-release/./primecount -t4 --test >/dev/null
        end=$(date +%s.%N)
        seconds_old=$(echo "scale=3; ($end - $start) / 1" | bc -l)
        echo "Seconds old code: $seconds_old"
        sleep 1

        start=$(date +%s.%N)
        build-curr-release/./primecount -t4 --test >/dev/null
        end=$(date +%s.%N)
        seconds_new=$(echo "scale=3; ($end - $start) / 1" | bc -l)
        echo "Seconds new code: $seconds_new"
        sleep 1

        limit=$(echo "scale=3; $seconds_old * $factor" | bc -l)
        new_code_is_fast=$(echo $seconds_new'<='$limit | bc -l)
        new_code_percent=$(echo "scale=1; 100 * $seconds_new / $seconds_old" | bc -l)

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
    factor=1.03

    # Test failure must be observed 3 times,
    # we try to avoid false negatives.
    for j in {1..3}
    do
        echo ""
        echo "=== Benchmark PrimePi(1e17) ===="
        echo ""

        total_seconds_new=0
        total_seconds_old=0

        echo "The test order is:"
        echo "1st current release - 2nd prev release."
        echo ""

        for i in {1..7}
        do
            seconds_new=$(build-curr-release/./primecount 1e17 -t4 --time | grep Seconds | cut -d' ' -f2)
            echo "Seconds new code: $seconds_new"
            total_seconds_new=$(echo "scale=3; ($total_seconds_new + $seconds_new) / 1" | bc -l)
            sleep 1
            seconds_old=$(build-prev-release/./primecount 1e17 -t4 --time | grep Seconds | cut -d' ' -f2)
            echo "Seconds old code: $seconds_old"
            echo ""
            total_seconds_old=$(echo "scale=3; ($total_seconds_old + $seconds_old) / 1" | bc -l)
            sleep 1
        done

        echo "Now reverse the test order:"
        echo "1st prev release - 2nd current release."
        echo ""

        for i in {1..7}
        do
            seconds_old=$(build-prev-release/./primecount 1e17 -t4 --time | grep Seconds | cut -d' ' -f2)
            echo "Seconds old code: $seconds_old"
            total_seconds_old=$(echo "scale=3; ($total_seconds_old + $seconds_old) / 1" | bc -l)
            sleep 1
            seconds_new=$(build-curr-release/./primecount 1e17 -t4 --time | grep Seconds | cut -d' ' -f2)
            echo "Seconds new code: $seconds_new"
            echo ""
            total_seconds_new=$(echo "scale=3; ($total_seconds_new + $seconds_new) / 1" | bc -l)
            sleep 1
        done

        echo "Total seconds old code: $total_seconds_old"
        echo "Total seconds new code: $total_seconds_new"
        echo ""

        limit=$(echo "scale=3; $total_seconds_old * $factor" | bc -l)
        new_code_is_fast=$(echo $total_seconds_new'<='$limit | bc -l)
        new_code_percent=$(echo "scale=1; 100 * $total_seconds_new / $total_seconds_old" | bc -l)

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
echo ""

echo "=== Old code (previous release) ==="
build-prev-release/./primecount --version
echo ""

echo "=== New code ==="
build-curr-release/./primecount --version

# Execute benchmark functions
benchmark_test_option
benchmark_pi_1e17
