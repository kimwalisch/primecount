# primecount build instructions

## Prerequisites

You need to have installed a C++ compiler which supports C++11 (or later) and CMake ≥ 3.4.

<table>
    <tr>
        <td><b>Debian/Ubuntu:</b></td>
        <td><code>sudo apt install g++ cmake</code></td>
    </tr>
    <tr>
        <td><b>Fedora:</b></td>
        <td><code>sudo dnf install gcc-c++ cmake</code></td>
    </tr>
    <tr>
        <td><b>openSUSE:</b></td>
        <td><code>sudo zypper install gcc-c++ cmake</code></td>
    </tr>
    <tr>
        <td><b>Arch Linux:</b></td>
        <td><code>sudo pacman -S gcc cmake</code></td>
    </tr>
</table>

## Unix-like OSes

Open a terminal, cd into the primecount directory and run:

```bash
cmake .
make -j
sudo make install
```

## macOS

On macOS the default C++ compiler that can be installed using ```xcode-select --install```
does not support OpenMP multi-threading. Hence I suggest installing an alternative
C++ compiler that supports OpenMP.

```bash
# Install C++ compiler with OpenMP support
brew install cmake llvm libomp

# Build primecount with OpenMP
LIBRARY_PATH=$(brew --prefix llvm)/lib CXX=$(brew --prefix llvm)/bin/clang++ cmake .
make -j
sudo make install
```

## MinGW/MSYS2 (Windows)

Open a terminal, cd into the primecount directory and run:

```bash
cmake -G "Unix Makefiles" .
make -j
sudo make install
```

## Microsoft Visual C++

First install [Visual Studio](https://visualstudio.microsoft.com/downloads/)
(includes CMake) on your Windows PC. Then go to the start menu, select Visual
Studio and open a **x64 Command Prompt**. Now cd into the primecount directory
and run the commands below:

```bash
# Use 'cmake -G' to find your Visual Studio version
cmake -G "Visual Studio 15 2017 Win64" .
cmake --build . --config Release

# Optionally install using Admin shell
cmake --build . --config Release --target install
```

## Run the tests

Open a terminal, cd into the primecount directory and run:

```bash
cmake -DBUILD_TESTS=ON .
make -j
ctest
```

## CMake configure options

By default the primecount binary, the static libprimecount and
libprimesieve will be built. The build options can be modified at
the configure step using e.g. ```cmake . -DBUILD_TESTS=ON```.

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
option(WITH_MPI            "Enable MPI support"                    OFF)
```

## Packaging primecount

When packaging primecount for e.g. a Linux distro it is best to change
a few of the default options.

primecount includes the libprimesieve dependency in its source tree and
libprimesieve is built by default. When packaging primecount it is better
to install libprimesieve using the package manager and not build
libprimesieve from source. You can achieve this using:

* ```cmake . -DBUILD_LIBPRIMESIEVE=OFF```

By default primecount builds the primecount binary and the static
libprimecount. Usually Linux distros don't want to package static
libraries. Hence you can build the primecount binary and the shared
libprimecount using:

* ```cmake . -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=OFF```

For performance reasons primecount uses the POPCNT instruction on all
CPU architectures that support it. On x86/x64 the POPCNT instruction was
added to Intel's and AMD's CPUs alongside the SSE4 instruction set in 2008.
If you need to support older x86/x64 CPUs you can disable POPCNT but this
will deteriorate performance by about 40%. Ideally POPCNT support should
only be disabled for x86/x64 since other CPU architectures
(e.g. ARM64, PPC64LE) do not have the issue that there exist CPUs with and
without the POPCNT instruction for the same instruction set.

* ```cmake . -DWITH_POPCNT=OFF```

## Man page regeneration

primecount includes an up to date man page at ```doc/primecount.1```.
That man page has been generated from ```doc/primecount.txt``` using
the ```a2x``` program from the ```asciidoc``` package. Usually when packaging
primecount it is recommended to regenerate the man page. In order to
regenerate the man page you need to install the ```asciidoc``` package and
then build primecount using ```cmake . -DBUILD_MANPAGE=ON```.

<table>
    <tr>
        <td><b>Debian/Ubuntu:</b></td>
        <td><code>sudo apt install asciidoc-base</code></td>
    </tr>
    <tr>
        <td><b>Fedora:</b></td>
        <td><code>sudo dnf install asciidoc</code></td>
    </tr>
    <tr>
        <td><b>openSUSE:</b></td>
        <td><code>sudo zypper install asciidoc</code></td>
    </tr>
    <tr>
        <td><b>Arch Linux:</b></td>
        <td><code>sudo pacman -S asciidoc</code></td>
    </tr>
</table>
