# libprimecount

primecount can be built as a static and shared C/C++ library for use in
other math projects. libprimecount has both a C API (```<primecount.h>``` header)
and a C++ API (```<primecount.hpp>``` header) so you are free to pick the one
that best fits your needs. The C API has been added to make it easier to write
bindings for other programming languages for libprimecount.

primecount's prime counting function implementation and nth prime function are
currently (March 2018) orders of magnitude faster than e.g. Mathematica, PARI/GP,
SageMath and SymPy. libprimecount is also very portable, it has been tested
successfully on a wide range of operating systems, compilers (GCC, Clang, MSVC)
and CPU architectures (x86, x64, ARM, ARM64, PowerPC, PP64, Sparc).

## Build instructions

You need to have installed a C++ compiler, cmake and make.

```sh
cmake . -DBUILD_SHARED_LIBS=ON
make -j
sudo make install
```

#### Run the tests

```sh
cmake . -DBUILD_TESTS=ON
make -j
make test
```

#### Maximum portability

By default primecount uses the ```POPCNT``` instruction in order to achieve the
best performance. As a drawback primecount won't work on CPUs that do not
have the ```POPCNT``` instruction e.g. all x86 CPUs built before 2010 do not
have the ```POPCNT``` instruction. If you require primecount to run on all CPUs
you have to disable ```POPCNT```:

```
cmake . -DWITH_POPCNT=OFF
```

#### CMake build options

Here are all available cmake configuration options:

```CMake
option(BUILD_PRIMECOUNT    "Build the primecount binary"           ON)
option(BUILD_LIBPRIMESIEVE "Build libprimesieve"                   ON)
option(BUILD_SHARED_LIBS   "Build the shared libprimecount"        OFF)
option(BUILD_STATIC_LIBS   "Build the static libprimecount"        ON)
option(BUILD_MANPAGE       "Regenerate man page using a2x program" OFF)
option(BUILD_TESTS         "Build the test programs"               OFF)

option(WITH_POPCNT         "Use the POPCNT instruction"            ON)
option(WITH_OPENMP         "Enable OpenMP multi-threading"         ON)
option(WITH_LIBDIVIDE      "Use libdivide.h"                       ON)
option(WITH_FLOAT128       "Use __float128 (requires libquadmath)" OFF)
option(WITH_MPI            "Enable MPI support"                    OFF)
```

## Linking

```sh
cc -O2 primes.c -lprimecount
c++ -O2 primes.cpp -lprimecount
```
