#!/bin/bash

# Usage: scripts/build_mingw64_arm64.sh
# Builds a primecount release binary that is statically linked
# and ready for distribution.

# Prerequisites:
# 1) Install MSYS2 x64 (or arm64 if available)
# 2) Open C:/msys64/clangarm64.exe
# 3) pacman -Syu (exit then run it again)
# 4) pacman -S mingw-w64-clang-aarch64-clang mingw-w64-clang-aarch64-openmp make cmake git unzip
# 5) cd primecount
# 6) scripts/build_mingw64_arm64.sh

# Exit if any error occurs
set -e

rm -rf build*

####################################################################

FULL_DATE=$(date +'%B %d, %Y')
YEAR=$(date +'%Y')

cd include
VERSION=$(grep "PRIMECOUNT_VERSION " primecount.hpp | cut -f2 -d'"')
cd ..

####################################################################

handle_error() {
    echo ""
    echo "Error: $1"
    exit 1
}

####################################################################

# The repo must no have any uncommited changes as we
# switch to another branch during the script.
git diff --exit-code > /dev/null || handle_error "repo must not have any uncommitted changes"

# Build primecount binary ##########################################

git pull
mkdir build-release-arm64
cd build-release-arm64

rm ../src/deleglise-rivat/S2_easy.cpp
rm ../src/gourdon/AC.cpp
clang++ -flto -static -O3 -DNDEBUG -D_WIN32_WINNT=0x0A00 -Wall -Wextra -pedantic -fopenmp -I ../include -I ../lib/primesieve/include ../lib/primesieve/src/*.cpp ../src/*.cpp ../src/lmo/*.cpp ../src/deleglise-rivat/*.cpp ../src/gourdon/*.cpp ../src/app/*.cpp -o primecount.exe -lPsapi
git checkout ..

strip primecount.exe

# Create a release zip archive
wget https://github.com/kimwalisch/primecount/releases/download/v6.1/primecount-6.1-win64.zip
unzip primecount-6.1-win64.zip -d primecount-$VERSION-win-arm64
rm primecount-6.1-win64.zip

echo ""
echo ""
echo "Old file size: $(ls -l --block-size=K primecount-$VERSION-win-arm64/primecount.exe)"
echo "New file size: $(ls -l --block-size=K primecount.exe)"
echo ""
echo ""

mv -f primecount.exe primecount-$VERSION-win-arm64
cd primecount-$VERSION-win-arm64
sed -i "1 s/.*/primecount $VERSION/" README.txt
sed -i "2 s/.*/$FULL_DATE/" README.txt
sed -i "3 s/.*/Copyright \(c\) 2013 - $YEAR, Kim Walisch\./" COPYING

# Verify sed has worked correctly
[ "$(sed -n '1p' < README.txt)" = "primecount $VERSION" ] || handle_error "failed updating README.txt"
[ "$(sed -n '2p' < README.txt)" = "$FULL_DATE" ] || handle_error "failed updating README.txt"
[ "$(sed -n '3p' < COPYING)" = "Copyright (c) 2013 - $YEAR, Kim Walisch." ] || handle_error "failed updating COPYING"

zip primecount-$VERSION-win-arm64.zip primecount.exe README.txt COPYING
cp primecount-$VERSION-win-arm64.zip ..

./primecount --test
echo ""
echo ""
./primecount 1e18 -s

cd ..

####################################################################

echo ""
echo "Release binary built successfully!"
