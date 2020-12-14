#!/bin/bash

# Usage: scripts/build_mingw64.sh
# Builds primecount and primecount-backup release binaries that
# are statically linked and ready for distribution.

# Exist if any error occurs
set -e

rm -rf build*

####################################################################

git checkout master
mkdir build-master
cd build-master
cmake .. -G "Unix Makefiles" -DCMAKE_CXX_FLAGS="-static -static-libgcc -static-libstdc++ -Wall -Wextra -pedantic -D_WIN32_WINNT=0x601" -DWITH_FLOAT128=ON
make -j8
rm primecount.exe

sed -i 's/-lkernel32.*/-lkernel32/g' CMakeFiles/primecount.dir/linklibs.rsp
sed -i 's/libgomp\.dll\.a/libgomp\.a/g' CMakeFiles/primecount.dir/linklibs.rsp

make
strip primecount.exe

####################################################################

cd ..
git checkout backup3
mkdir build-backup3
cd build-backup3
cmake .. -G "Unix Makefiles" -DCMAKE_CXX_FLAGS="-static -static-libgcc -static-libstdc++ -Wall -Wextra -pedantic -D_WIN32_WINNT=0x601" -DWITH_FLOAT128=ON
make -j8
rm primecount.exe

sed -i 's/-lkernel32.*/-lkernel32/g' CMakeFiles/primecount.dir/linklibs.rsp
sed -i 's/libgomp\.dll\.a/libgomp\.a/g' CMakeFiles/primecount.dir/linklibs.rsp

make
strip primecount.exe
