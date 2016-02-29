#!/bin/sh

# Usage: ./build.sh
# Script which automates building primecount and libprimecount.
#
# What this script does:
#
# 1) Download primesieve library
# 2) Build primesieve library using: ./configure && make
# 3) Build primecount library using: ./configure && make
#
# Lots of hacks are needed because:
#
# 1) We build primecount without first installing libprimesieve.
# 2) We want to build a static primecount binary.

CONFIGURE_OPTIONS="$1"

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
    cd primesieve-*
    ./configure
    make -j8
    # On Windows we must install libprimesive
    if [ "$(uname | egrep -i 'windows|cygwin|mingw|msys')" != "" ]
    then
        make install
        export LIBRARY_PATH=/usr/local/lib:$LIBRARY_PATH
    fi
    cd ..
fi

# Generate ./configure script, requires Autotools
if [ ! -f ./configure ]
then
    ./autogen.sh
fi

# Patch ./configure script to continue even
# if libprimesieve is not installed
if [ "$(grep 'libprimesieve is missing' configure)" != "" ]
then
    sed '/libprimesieve is missing/c\
    true;
    ' configure > configure.tmp
    mv -f configure.tmp configure
    chmod +x configure
fi

# Generate Makefile using ./configure
if [ ! -f ./Makefile ]
then
    ./configure $CONFIGURE_OPTIONS CXXFLAGS="-O2 -I$(echo primesieve-*/include)"
fi

# Patch Makefile to build primecount binary which links statically
# against libprimecount and libprimesieve
if [ "$(grep libprimesieve.a Makefile)" = "" ]
then
    sed '/primecount_DEPENDENCIES = libprimecount.la/c\
    primecount_DEPENDENCIES = libprimecount.la primesieve-*/.libs/libprimesieve.a
    ' Makefile > Makefile.tmp
    mv -f Makefile.tmp Makefile

    sed '/primecount_LDADD = libprimecount.la/c\
    primecount_LDADD = libprimecount.la primesieve-*/.libs/libprimesieve.a
    ' Makefile > Makefile.tmp
    mv -f Makefile.tmp Makefile

    chmod +x Makefile
fi

# Build both static and shared libprimecount
make libprimecount.la -j8

# Build statically linked primecount binary
make primecount$(grep 'EXEEXT =' Makefile | cut -f3 -d' ') LDFLAGS="-static" -j8
