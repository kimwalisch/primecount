#!/bin/bash

# Usage: scripts/build_mingw64_x64.sh
# Builds a primecount release binary that is statically linked
# and ready for distribution.

# === Prerequisites x64 ===
# 1) Install MSYS2 x64
# 2) pacman -Syu (exit then run it again)
# 3) pacman -S --needed base-devel mingw-w64-x86_64-toolchain mingw-w64-x86_64-cmake zip unzip git
# 4) git clone https://github.com/kimwalisch/primecount.git
# 5) scripts/build_mingw64_x64.sh

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
mkdir build-release
cd build-release
cmake .. -G "Unix Makefiles" -DCMAKE_CXX_FLAGS="-flto -static -static-libgcc -static-libstdc++ -Wall -Wextra -pedantic -D_WIN32_WINNT=0x601" -DWITH_FLOAT128=ON
make -j8
rm primecount.exe

# Remove unnecessary libraries for linking,
# keep only GCC libraries + kernel32.
sed -i 's/-lkernel32.*/-lkernel32/g' CMakeFiles/primecount.dir/linklibs.rsp
sed -i 's/libgomp\.dll\.a/libgomp\.a/g' CMakeFiles/primecount.dir/linklibs.rsp

# Verify that sed has worked correctly,
# last word should be -lkernel32.
[ "$(grep -o '[^ ]\+$' CMakeFiles/primecount.dir/linklibs.rsp)" = "-lkernel32" ] || handle_error "failed updating linklibs.rsp"

make
strip primecount.exe

# Create a release zip archive
wget https://github.com/kimwalisch/primecount/releases/download/v7.11/primecount-7.11-win-x64.zip
unzip primecount-7.11-win-x64.zip -d primecount-$VERSION-win-x64
rm primecount-7.11-win-x64.zip

echo ""
echo ""
echo "Old file size: $(ls -l --block-size=K primecount-$VERSION-win-x64/primecount.exe)"
echo "New file size: $(ls -l --block-size=K primecount.exe)"
echo ""
echo ""

mv -f primecount.exe primecount-$VERSION-win-x64
cd primecount-$VERSION-win-x64
sed -i "1 s/.*/primecount $VERSION/" README.txt
sed -i "2 s/.*/$FULL_DATE/" README.txt
sed -i "3 s/.*/Copyright \(c\) 2013 - $YEAR, Kim Walisch\./" COPYING

# Verify sed has worked correctly
[ "$(sed -n '1p' < README.txt)" = "primecount $VERSION" ] || handle_error "failed updating README.txt"
[ "$(sed -n '2p' < README.txt)" = "$FULL_DATE" ] || handle_error "failed updating README.txt"
[ "$(sed -n '3p' < COPYING)" = "Copyright (c) 2013 - $YEAR, Kim Walisch." ] || handle_error "failed updating COPYING"

zip primecount-$VERSION-win-x64.zip primecount.exe README.txt COPYING
cp primecount-$VERSION-win-x64.zip ..

./primecount --test
echo ""
echo ""
./primecount 1e18 -s

cd ..

####################################################################

echo ""
echo "Release binary built successfully!"
