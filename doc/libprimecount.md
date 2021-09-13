# libprimecount

primecount can be built as a static and shared C/C++ library for use in other
math projects. libprimecount has both a C API (```<primecount.h>``` header) and
a C++ API (```<primecount.hpp>``` header) so you are free to pick the one that
best fits your needs. The C API has been added to make it easier to write
libprimecount bindings for other programming languages.

primecount's prime counting function implementation and nth prime function are
currently (March 2018) orders of magnitude faster than e.g. Mathematica, PARI/GP,
SageMath and SymPy. libprimecount is also very portable, it has been tested
successfully on a wide range of operating systems, compilers (GCC, Clang, MSVC)
and CPU architectures (x86, x64, ARM, ARM64, PowerPC, PP64, Sparc).

# C example

The C example below counts the primes ≤ 1000 and prints the result to the screen.
Note that primecount is multi-threaded by default, it uses all available CPU
cores if the input is sufficiently large.

```C
#include <primecount.h>
#include <stdio.h>

int main()
{
    int64_t pix = primecount_pi(1000);
    printf("primes below 1000 = %ld\n", pix);

    return 0;
}
```

Compile using:

```sh
cc -O2 primes.c -lprimecount
```

# C API reference

Include the ```<primecount.h>``` header to use primecount's C API.
All functions that are part of primecount's C API return ```-1``` in case an
error occurs and print the corresponding error message to the standard error
stream.

```C
// Count the number of primes <= x
int64_t primecount_pi(int64_t x);

// Count the number of primes <= x (supports 128-bit)
int primecount_pi_str(const char* x, char* res, size_t len);

// Find the nth prime e.g.: nth_prime(25) = 97
int64_t primecount_nth_prime(int64_t n);

// Count the numbers <= x that are not divisible by any of the first a primes
int64_t primecount_phi(int64_t x, int64_t a);
```

Please see [primecount.h](https://github.com/kimwalisch/primecount/blob/master/include/primecount.h)
for more information.

# C++ example


The C++ example below counts the primes ≤ 1000 and prints the result to the screen.
Note that primecount is multi-threaded by default, it uses all available CPU
cores if the input is sufficiently large.

```C++
#include <primecount.hpp>
#include <iostream>

int main()
{
    int64_t pix = primecount::pi(1000);
    std::cout << "primes below 1000 = " << pix << std::endl;

    return 0;
}
```

Compile using:

```sh
c++ -O2 primes.cpp -lprimecount
```

# C++ API reference

Include the ```<primecount.hpp>``` header to use primecount's C++ API.
All functions that are part of primecount's C++ API throw a
```primecount_error``` exception (which is derived from
```std::exception```) in case an error occurs.

```C++
// Count the number of primes <= x
int64_t primecount::pi(int64_t x);

// Count the number of primes <= x (supports 128-bit)
std::string primecount::pi(const std::string& x);

// Find the nth prime e.g.: nth_prime(25) = 97
int64_t primecount::nth_prime(int64_t n);

// Count the numbers <= x that are not divisible by any of the first a primes
int64_t primecount::phi(int64_t x, int64_t a);
```

Please see [primecount.hpp](https://github.com/kimwalisch/primecount/blob/master/include/primecount.hpp)
for more information.

# Build instructions

You need to have [installed](BUILD.md#prerequisites) a C++ compiler, cmake and make.

```sh
cmake . -DBUILD_SHARED_LIBS=ON
make -j
sudo make install
```

# Run the tests

```sh
cmake . -DBUILD_TESTS=ON
make -j
make test
```

# Maximum portability

By default libprimecount uses the ```POPCNT``` instruction in order to achieve the
best performance. As a drawback libprimecount won't work on CPUs that do not
have the ```POPCNT``` instruction e.g. all x86 CPUs built before 2010 do not
have the ```POPCNT``` instruction. If you require libprimecount to run on all CPUs
you have to disable ```POPCNT```:

```
cmake . -DWITH_POPCNT=OFF
```

# CMake build options

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
option(WITH_DIV32          "Use 32-bit division instead of 64-bit division whenever possible" ON)
option(WITH_FLOAT128       "Use __float128 (requires libquadmath)" OFF)
option(WITH_JEMALLOC       "Use jemalloc allocator"                OFF)
```
