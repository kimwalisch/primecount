#!/bin/bash

# Usage: scripts/build_mingw64_arm64.sh
# Builds a primecount release binary that is statically linked
# and ready for distribution.

# === Prerequisites arm64 ===
# 1) Install a trial version of both Parallels & Windows on a MacBook ARM64.
# 2) No need to purchase/register Parallels & Windows, keep using the trial version.
# 3) Install MSYS2 x64 (or arm64 if available)
# 4) Open C:/msys64/clangarm64.exe
# 5) pacman -Syu (exit then run it again)
# 6) pacman -S mingw-w64-clang-aarch64-clang mingw-w64-clang-aarch64-openmp make git zip unzip
# 7) git clone https://github.com/kimwalisch/primecount.git
# 8) scripts/build_mingw64_arm64.sh

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
rm ../src/deleglise-rivat/S2_hard_multiarch_arm_sve.cpp
rm ../src/deleglise-rivat/S2_hard_multiarch_avx512.cpp
rm ../src/gourdon/AC.cpp
rm ../src/gourdon/D_multiarch_arm_sve.cpp
rm ../src/gourdon/D_multiarch_avx512.cpp

mkdir build_primesieve
cd build_primesieve
clang++ -c -I../../lib/primesieve/include -I../../lib/primesieve/src \
  -O3 -flto -static -Wall -Wextra -pedantic \
  -DENABLE_MULTIARCH_ARM_SVE -DNDEBUG -D_WIN32_WINNT=0x0A00 \
  ../../lib/primesieve/src/*.cpp ../../lib/primesieve/src/arch/arm/sve.cpp

cd ..
mkdir build_primecount
cd build_primecount
clang++ -c -I../../include -I../../src -I../../lib/primesieve/include \
  -O3 -flto -fopenmp -static -Wall -Wextra -pedantic \
  -DENABLE_MULTIARCH_ARM_SVE -DNDEBUG -D_WIN32_WINNT=0x0A00 \
  ../../src/*.cpp ../../src/lmo/*.cpp ../../src/deleglise-rivat/*.cpp \
  ../../src/gourdon/*.cpp ../../src/arch/arm/sve.cpp ../../src/app/*.cpp

cd ..
clang++ -O3 -flto -fopenmp -static -Wall -Wextra -pedantic -DENABLE_MULTIARCH_ARM_SVE -DNDEBUG -D_WIN32_WINNT=0x0A00 \
  build_primesieve/*.o build_primecount/*.o -o primecount -lPsapi

strip primecount.exe
git checkout ..

# Create a release zip archive
wget https://github.com/kimwalisch/primecount/releases/download/v7.11/primecount-7.11-win-arm64.zip
unzip primecount-7.11-win-arm64.zip -d primecount-$VERSION-win-arm64
rm primecount-7.11-win-arm64.zip

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
