#!/bin/sh

# Usage: ./build.sh
# Script which automates building primecount.
# Prerequisites: cmake & make.

BUILD_OPTIONS="$@"
CPU_CORES=$(nproc --all 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 8)

# Exit on any error
set -e

# Download primesieve
if [ ! -f ./primesieve-latest.tar.gz ]
then
    wget      https://dl.bintray.com/kimwalisch/primesieve/primesieve-latest.tar.gz || \
    curl -LkO https://dl.bintray.com/kimwalisch/primesieve/primesieve-latest.tar.gz
fi

# Build libprimesieve
if [ ! -f primesieve*/.libs/libprimesieve.a ]
then
    tar xvf primesieve-latest.tar.gz
    mv primesieve-*/ primesieve-latest
    cd primesieve-latest
    ./configure --disable-shared
    make -j$CPU_CORES
    cd ..
fi

# Generate Makefile
cmake "$BUILD_OPTIONS" .

# Build statically linked primecount binary
make -j$CPU_CORES
