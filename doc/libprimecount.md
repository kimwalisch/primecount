# libprimecount

primecount can be built as a static and shared C++ library for use in
other math projects. primecount's prime counting function implementation
and nth prime function are currently (March 2018) orders of magnitude
faster than e.g. Mathematica, PARI/GP, SageMath and SymPy.

libprimecount is also very portable, it has been tested successfully on
a wide range of operating systems, compilers (GCC, Clang, MSVC) and CPU
architectures (x86, x64, ARM, ARM64, PowerPC, PP64, Sparc).

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

## C++ API

All functions are multi-threaded by default.

```C++
#include <primecount.hpp>

/// Count the primes <= x
int64_t primecount::pi(int64_t x);

/// 128-bit prime counting function.
/// @param expr  Integer arithmetic expression e.g. "1000", "10^22"
/// @pre   expr  <= 10^31 on 64-bit systems
///        expr    < 2^63 on 32-bit systems
std::string primecount::pi(const std::string& expr);

/// Find the nth prime using a combination of the prime
/// counting function and the sieve of Eratosthenes.
/// Run time: O(x^(2/3) / (log x)^2)
/// Memory usage: O(x^(1/2))
///
int64_t nth_prime(int64_t n);

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
///
int64_t phi(int64_t x, int64_t a);
```

## Usage example

Below is an example program that counts the primes below 1000.

```C++
#include <primecount.hpp>
#include <iostream>

int main()
{
    int64_t primes = primecount::pi(1000);
    std::cout << "primes below 1000 = " << primes << std::endl;
  
    return 0;
}
```

## Linking

```sh
c++ -O2 primes.cpp -lprimecount
```
