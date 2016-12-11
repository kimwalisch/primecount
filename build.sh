#!/bin/sh

# Usage: ./build.sh
# Script which automates building primecount.
# Prerequisites: make & cmake.

CMAKE_OPTIONS="$@"
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

# MinGW/MSYS (Windows) requires special cmake command
uname 2>/dev/null | grep -i MinGW >/dev/null
if [ $? -eq 0 ];
then
    cmake -G "MSYS Makefiles" "$CMAKE_OPTIONS" .
    make -j$CPU_CORES
    exit 0
fi

# Generate Makefile
cmake "$CMAKE_OPTIONS" .

# Build primecount
make -j$CPU_CORES
