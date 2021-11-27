# libprimecount

primecount can be built as a static and shared C/C++ library for use in other
math projects. libprimecount has both a C API (```<primecount.h>``` header) and
a C++ API (```<primecount.hpp>``` header) so you are free to pick the one that
best fits your needs. The C API has been added to make it easier to write
libprimecount bindings for other programming languages.

primecount's prime counting function implementation and nth prime function are
generally orders of magnitude faster than other publicly available prime counting
function implementations. As of 2021 libprimecount is the only prime counting
function library that I am aware of that supports multi-threading. libprimecount
is also very portable, it has been tested successfully on a wide range of
operating systems, compilers (GCC, Clang, MSVC) and CPU architectures (x86, x64,
ARM, ARM64, PowerPC, PP64, Sparc).

libprimecount has been integrated into a few other computer algebra systems:
Since 2021
[Mathematica's PrimePi(x)](https://reference.wolfram.com/language/ref/PrimePi.html)
uses libprimecount under the hood, primecount is available as an
optional package in
[SageMath](https://doc.sagemath.org/html/en/reference/spkg/primecount.html) and
there are bindings available for a few
[other programming languages](https://github.com/kimwalisch/primecount#bindings-for-other-languages).
If your company uses primecount then I would greatly appreciate if you would
[sponsor](https://github.com/sponsors/kimwalisch) its maintenance (or
[donate](https://github.com/sponsors/kimwalisch?frequency=one-time&sponsor=kimwalisch)).
The development of primecount has personally cost me a lot of money, as I
frequently needed to run extensive benchmarks on a wide variety of high end
cloud servers.

# Installation

* [Install libprimecount using package manager](https://github.com/kimwalisch/primecount#installation)
* [Build libprimecount from source](#build-instructions)

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
cc -O3 primes.c -o primes -lprimecount
```

If you have [built libprimecount yourself](#Build-instructions),
then the default installation path is usually ```/usr/local/lib```. Running
the ```ldconfig``` program after ```make install``` ensures that Linux's dynamic
linker/loader will find the shared primecount library when you execute your program.
However, some OSes are missing the ```ldconfig``` program or ```ldconfig``` does
not include ```/usr/local/lib``` by default. In these cases you need to export
some environment variables:

```sh
export LIBRARY_PATH=/usr/local/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=/usr/local/include:$C_INCLUDE_PATH
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
c++ -O3 primes.cpp -o primes -lprimecount
```

If you have [built libprimecount yourself](#Build-instructions),
then the default installation path is usually ```/usr/local/lib```. Running
the ```ldconfig``` program after ```make install``` ensures that Linux's dynamic
linker/loader will find the shared primecount library when you execute your program.
However, some OSes are missing the ```ldconfig``` program or ```ldconfig``` does
not include ```/usr/local/lib``` by default. In these cases you need to export
some environment variables:

```sh
export LIBRARY_PATH=/usr/local/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
export CPLUS_INCLUDE_PATH=/usr/local/include:$CPLUS_INCLUDE_PATH
```

# Build instructions

You need to have [installed a C++ compiler, cmake and make](BUILD.md#prerequisites). By default,
only the static libprimecount is built, but for doing development with libprimecount you will
most likely want to use the shared libprimecount. Hence we use the ```-DBUILD_SHARED_LIBS=ON```
option to enable building libprimecount as a shared library.

```sh
cmake . -DBUILD_SHARED_LIBS=ON
make -j
sudo make install
sudo ldconfig
```

* [Detailed build instructions](BUILD.md#primecount-build-instructions)

# Maximum portability

For performance reasons primecount uses the ```POPCNT``` instruction on all CPU architectures that
support it. On x86/x64 the ```POPCNT``` instruction was added to Intel's and AMD's CPUs alongside the
SSE4 instruction set in 2008. If you need to support older x86/x64 CPUs you can disable ```POPCNT``` but
this will deteriorate performance by about 50%. Note that disabling ```POPCNT``` only has an effect on
x86/x64, on other CPU architectures ```POPCNT``` is always used if it is available (as this generally
does not cause any issues).

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

option(WITH_POPCNT          "Use the POPCNT instruction"            ON)
option(WITH_LIBDIVIDE       "Use libdivide.h"                       ON)
option(WITH_OPENMP          "Enable OpenMP multi-threading"         ON)
option(WITH_DIV32           "Use 32-bit division instead of 64-bit division whenever possible" ON)
option(WITH_MSVC_CRT_STATIC "Link primecount.lib with /MT instead of the default /MD" OFF)
option(WITH_FLOAT128        "Use __float128 (requires libquadmath)" OFF)
option(WITH_JEMALLOC        "Use jemalloc allocator"                OFF)
```
