:: libatomic on Windows does not yet support 128-bit integers,
:: hence currently the build fails with libatomic linker errors.

:: After a new LLVM/Clang installation I should update the file
:: C:\Program Files\LLVM\lib\clang\18\include\immintrin.h and unconditionally
:: include all MMX, SSE, POPCNT, BMI, BMI2, AVX, AVX and AVX512 headers.

del /Q ..\src\deleglise-rivat\S2_easy.cpp
del /Q ..\src\deleglise-rivat\S2_hard_multiarch_arm_sve.cpp
del /Q ..\src\gourdon\AC.cpp
del /Q ..\src\gourdon\D_multiarch_arm_sve.cpp
clang++ -I../include -I../lib/primesieve/include -O3 -mpopcnt -fopenmp -Wall -Wextra -pedantic -DNDEBUG -DENABLE_MULTIARCH_AVX512_BW -DENABLE_MULTIARCH_AVX512_VBMI2 -DENABLE_MULTIARCH_AVX512_VPOPCNT ../lib/primesieve/src/*.cpp ../lib/primesieve/src/arch/x86/*.cpp ../src/*.cpp ../src/arch/x86/*.cpp ../src/lmo/*.cpp ../src/deleglise-rivat/*.cpp ../src/gourdon/*.cpp ../src/app/*.cpp -o primecount.exe "C:\Program Files\LLVM\lib\clang\18\lib\windows\clang_rt.builtins-x86_64.lib"
git checkout ..\src\deleglise-rivat
git checkout ..\src\gourdon
