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

# Download and build primesieve
wget https://github.com/kimwalisch/primesieve/archive/master.tar.gz || curl -OLk https://github.com/kimwalisch/primesieve/archive/master.tar.gz
tar xvf master.tar.gz
cd primesieve-master
./autogen.sh
./configure --prefix=$(pwd)/..
make -j8
make install
cd ..

# Patch Makefile.am for static linking libprimesieve
sed 's#primecount_LDADD = libprimecount.la#primecount_LDADD = libprimecount.la lib/libprimesieve.a#g' Makefile.am > Makefile.tmp
mv -f Makefile.tmp Makefile.am

# Build primecount-mpi
./autogen.sh
./configure --disable-shared LDFLAGS=-L$(pwd)/lib
make -j8
