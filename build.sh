#!/bin/sh

# Usage: ./build.sh
# Builds a static primecount binary, automatically downloads and
# builds the primesieve library (dependency).

# Exit on any error
set -e

CONFIGURE_OPTIONS="$1"

# Download primesieve
if [ ! -f ./primesieve-latest.tar.gz ]
then
    wget      http://dl.bintray.com/kimwalisch/primesieve/primesieve-latest.tar.gz || \
    curl -LkO http://dl.bintray.com/kimwalisch/primesieve/primesieve-latest.tar.gz
fi

# Build libprimesieve
if [ ! -f lib/libprimesieve.a ]
then
    tar xvf primesieve-latest.tar.gz
    cd primesieve-*
    ./configure --prefix=$(pwd)/..
    make -j8
    make install
    cd ..
fi

# Generate configure script, requires GNU Autotools
if [ ! -f ./configure ]
then
    ./autogen.sh
fi

# configure primecount-mpi
if [ ! -f ./Makefile ]
then
    ./configure $CONFIGURE_OPTIONS LDFLAGS="-static -L$(pwd)/lib"
fi

# Build primecount-mpi
make -j8
