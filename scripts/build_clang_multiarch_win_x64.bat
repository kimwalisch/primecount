:: libatomic on Windows does not yet support 128-bit integers, hence for the
:: time being we need to use our -DENABLE_INT128_OPENMP_PATCH workaround
:: which forces OpenMP to use critical sections instead of atomics for
:: 128-bit integers. This workaround is slightly less performant than atomics.
:: Hence the workaruond should be removed once libatomic supports 128-bit
:: integers.

:: After a new LLVM/Clang installation I should update the file
:: C:\Program Files\LLVM\lib\clang\18\include\immintrin.h and unconditionally
:: include all MMX, SSE, POPCNT, BMI, BMI2, AVX, AVX and AVX512 headers.

del /Q ..\src\deleglise-rivat\S2_easy.cpp ..\src\gourdon\AC.cpp
clang++ -I../include -I../lib/primesieve/include -O3 -mpopcnt -fopenmp -Wall -Wextra -pedantic -DNDEBUG -DENABLE_INT128_OPENMP_PATCH -DENABLE_MULTIARCH_AVX512_VBMI2 -DENABLE_MULTIARCH_AVX512_VPOPCNT ../lib/primesieve/src/*.cpp ../lib/primesieve/src/x86/*.cpp ../src/*.cpp ../src/x86/*.cpp ../src/lmo/*.cpp ../src/deleglise-rivat/*.cpp ../src/gourdon/*.cpp ../src/app/*.cpp -o primecount.exe "C:\Program Files\LLVM\lib\clang\18\lib\windows\clang_rt.builtins-x86_64.lib"
git checkout ..\src\deleglise-rivat
git checkout ..\src\gourdon
