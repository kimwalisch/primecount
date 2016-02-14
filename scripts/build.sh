#!/bin/sh

# Usage: ./build.sh
# Builds a static primecount binary, automatically downloads and
# builds the primesieve library (dependency).

# Exit on any error
set -e

# Execute in base directory
if [ "$(basename $(pwd))" == "scripts" ]
then
    cd ..
fi

# Download primesieve
if [ ! -f ./primesieve-master.tar.gz ]
then
    wget -O primesieve-master.tar.gz https://github.com/kimwalisch/primesieve/archive/master.tar.gz || \
    curl -o primesieve-master.tar.gz -Lk https://github.com/kimwalisch/primesieve/archive/master.tar.gz
fi

# Build libprimesieve
if [ ! -f lib/libprimesieve.a ]
then
    tar xvf primesieve-master.tar.gz
    cd primesieve-master
    ./autogen.sh
    ./configure --prefix=$(pwd)/..
    make -j8
    make install
    cd ..
fi

# configure primecount-mpi
if [ ! -f Makefile ]
then
    if [ "$(grep libprimesieve.a Makefile.am)" = "" ]
    then
        # Patch Makefile.am for static linking libprimesieve
        sed 's#primecount_LDADD = libprimecount.la#primecount_LDADD = libprimecount.la lib/libprimesieve.a#g' Makefile.am > Makefile.tmp
        mv -f Makefile.tmp Makefile.am
    fi
    ./autogen.sh
    ./configure --disable-shared LDFLAGS=-L$(pwd)/lib
fi

# Build primecount-mpi
make -j8
