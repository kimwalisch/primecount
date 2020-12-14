#!/bin/bash

# Usage: scripts/build_mingw64.sh
# Builds primecount and primecount-backup release binaries that
# are statically linked and ready for distribution.

# Exist if any error occurs
set -e

rm -rf build*

####################################################################

FULL_DATE=$(date +'%B %d, %Y')
YEAR=$(date +'%Y')

cd include
VERSION=$(grep "PRIMECOUNT_VERSION " primecount.hpp | cut -f2 -d'"')
cd ..

# Build primecount binary ##########################################

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

# Create a release zip archive
wget https://github.com/kimwalisch/primecount/releases/download/v6.1/primecount-6.1-win64.zip
unzip primecount-6.1-win64.zip -d primecount-$VERSION-win-x64
rm primecount-6.1-win64.zip
mv -f primecount.exe primecount-$VERSION-win-x64

cd primecount-$VERSION-win-x64
sed -i "1 s/.*/primecount $VERSION/" README.txt
sed -i "2 s/.*/$FULL_DATE/" README.txt
sed -i "3 s/.*/Copyright \(c\) 2013 - $YEAR, Kim Walisch\./" COPYING

zip primecount-$VERSION-win-x64.zip primecount.exe README.txt COPYING
cp primecount-$VERSION-win-x64.zip ..
cd ..

# Build primecount-backup binary ###################################

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

# Create a release zip archive
wget https://github.com/kimwalisch/primecount/releases/download/v6.0-backup/primecount-backup-6.0-win64.zip
unzip primecount-backup-6.0-win64.zip -d primecount-backup-$VERSION-win-x64
rm primecount-backup-6.0-win64.zip

mv -f primecount.exe primecount-backup-$VERSION-win-x64
cd primecount-backup-$VERSION-win-x64
sed -i "1 s/.*/primecount-backup $VERSION/" README.txt
sed -i "2 s/.*/$FULL_DATE/" README.txt
sed -i "3 s/.*/Copyright \(c\) 2013 - $YEAR, Kim Walisch\./" COPYING

zip primecount-backup-$VERSION-win-x64.zip primecount.exe README.txt COPYING worktodo.sh worktodo.txt
cp primecount-backup-$VERSION-win-x64.zip ..
cd ..
