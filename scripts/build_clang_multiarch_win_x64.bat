REM After a new LLVM/Clang installation I should update the File
REM C:\Program Files\LLVM\lib\clang\18\include\immintrin.h and unconditionally
REM include all MMX, SSE, POPCNT, BMI, BMI2, AVX, AVX and AVX512 headers.

del /Q ..\src\deleglise-rivat\S2_easy.cpp ..\src\gourdon\AC.cpp
clang++ -I../include -I../lib/primesieve/include -O3 -mpopcnt -fopenmp -Wall -Wextra -pedantic -DNDEBUG -DENABLE_INT128_STL_PATCH -DENABLE_MULTIARCH_AVX512 -DENABLE_MULTIARCH_AVX512_BMI2 ../lib/primesieve/src/*.cpp ../src/*.cpp ../src/lmo/*.cpp ../src/deleglise-rivat/*.cpp ../src/gourdon/*.cpp ../src/app/*.cpp -o primecount.exe "C:\Program Files\LLVM\lib\clang\18\lib\windows\clang_rt.builtins-x86_64.lib"
git checkout ..
