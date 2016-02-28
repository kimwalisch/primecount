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
    wget      https://dl.bintray.com/kimwalisch/primesieve/primesieve-latest.tar.gz || \
    curl -LkO https://dl.bintray.com/kimwalisch/primesieve/primesieve-latest.tar.gz
fi

# Build libprimesieve
if [ ! -f lib/libprimesieve.a ]
then
    tar xvf primesieve-latest.tar.gz
    cd primesieve-*
    ./configure
    make -j8
    cd ..
fi

# Generate configure script, requires GNU Autotools
if [ ! -f ./configure ]
then
    # Patch configure.ac for static linking libprimesieve
    sed 's/AC_SEARCH_LIBS(\[primesieve/#AC_SEARCH_LIBS(\[primesieve/g' configure.ac > configure.tmp
    mv -f configure.tmp configure.ac
    ./autogen.sh
fi

# configure primecount-mpi
if [ ! -f ./Makefile ]
then
    if [ "$(grep libprimesieve.a Makefile.am)" = "" ]
    then
        # Patch Makefile.am for static linking libprimesieve
        sed 's#primecount_LDADD = libprimecount.la#primecount_LDADD = libprimecount.la primesieve*/.libs/libprimesieve.a#g' Makefile.am > Makefile.tmp
        mv -f Makefile.tmp Makefile.am
    fi
    ./configure $CONFIGURE_OPTIONS CXXFLAGS="-O2 -Iprimesieve*/include"
fi

make libprimecount.la -j8
make primecount LDFLAGS="-static" -j8
