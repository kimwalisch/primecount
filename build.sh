#!/bin/sh

# Usage: ./build.sh
# Script which automates building primecount.
# Prerequisites: make & cmake.

CMAKE_OPTIONS="$@"
CPU_CORES=$(nproc --all 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 8)
MINGW=$(uname 2>/dev/null | grep -i MinGW)

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
    ./configure
    make -j$CPU_CORES
    cd ..
fi

if [ ! -z "$MINGW" ]
then
    # MinGW/MSYS (Windows) requires special cmake command
    cmake -G "MSYS Makefiles" "$CMAKE_OPTIONS" .
    make -j$CPU_CORES
    exit 0
fi

# Generate Makefile
cmake "$CMAKE_OPTIONS" .

# Build primecount
make -j$CPU_CORES
