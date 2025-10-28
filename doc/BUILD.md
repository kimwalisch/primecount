# primecount build instructions

# Contents

* [Prerequisites](#prerequisites)
* [Linux & Unix-like OSes](#linux--unix-like-oses)
* [macOS](#macos)
* [MinGW/MSYS2 (Windows)](#mingwmsys2-windows)
* [Microsoft Visual C++](#microsoft-visual-c)
* [Emscripten/WebAssembly](#emscriptenwebassembly)
* [Run the tests](#run-the-tests)
* [CMake configure options](#cmake-configure-options)
* [Packaging primecount](#packaging-primecount)
* [Man page regeneration](#man-page-regeneration)

# Prerequisites

You need to have installed a C++ compiler which supports C++11 (or later) and CMake ≥ 3.4.

<table>
    <tr>
        <td><b>Arch Linux:</b></td>
        <td><code>sudo pacman -S gcc cmake</code></td>
    </tr>
    <tr>
        <td><b>Debian/Ubuntu:</b></td>
        <td><code>sudo apt install g++ cmake</code></td>
    </tr>
    <tr>
        <td><b>Fedora:</b></td>
        <td><code>sudo dnf install gcc-c++ cmake</code></td>
    </tr>
    <tr>
        <td><b>macOS:</b></td>
        <td><code>brew install cmake libomp</code></td>
    </tr>
    <tr>
        <td><b>openSUSE:</b></td>
        <td><code>sudo zypper install gcc-c++ cmake</code></td>
    </tr>
</table>

# Linux & Unix-like OSes

Open a terminal, cd into the primecount directory and run:

```bash
cmake .
cmake --build . --parallel
sudo cmake --install .
sudo ldconfig
```

# macOS

On macOS the default AppleClang C/C++ compiler can be installed using
```xcode-select --install```. If you have installed the [homebrew package
manager](https://brew.sh) then AppleClang is usually already installed on your Mac.

Open a terminal, cd into the primecount directory and run:

```bash
# Install CMake & the OpenMP library
brew install cmake libomp

CXXFLAGS="-I$(brew --prefix libomp)/include" LDFLAGS="-L$(brew --prefix libomp)/lib" cmake .
cmake --build . --parallel
sudo cmake --install .
```

# MinGW/MSYS2 (Windows)

Open a terminal, cd into the primecount directory and run:

```bash
cmake -G "Unix Makefiles" .
cmake --build . --parallel
```

# Microsoft Visual C++

First install [Visual Studio](https://visualstudio.microsoft.com/downloads/)
(includes CMake) on your Windows PC. Then go to the start menu, select Visual
Studio and open a **x64 Command Prompt**. Now cd into the primecount directory
and run the commands below:

```bash
# Use 'cmake -G' to find your Visual Studio version
cmake -G "Visual Studio 17 2022" .
cmake --build . --config Release

# Optionally install using Admin shell
cmake --install . --config Release
```

# Emscripten/WebAssembly

Using the Emscripten compiler you can compile the primecount C/C++ library to WebAssembly:

```bash
# Install the Emscripten compiler
git clone https://github.com/emscripten-core/emsdk.git
cd emsdk
./emsdk install latest
./emsdk activate latest
source emsdk_env.sh

# Compile primecount to WebAssembly
git clone https://github.com/kimwalisch/primecount.git
cd primecount
emcmake cmake .
emmake make -j4

# Run the primecount WebAssembly binary
node ./primecount.js 1e15 --status
```

# Run the tests

Open a terminal, cd into the primecount directory and run:

```bash
cmake . -DBUILD_TESTS=ON
cmake --build . --parallel
ctest
```

For developers hacking on primecount's source code the
[test/README.md](../test/README.md) document contains more information
about primecount testing such as testing in debug mode and testing
using GCC/Clang sanitizers.

# CMake configure options

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

option(WITH_LIBDIVIDE       "Use libdivide.h"                       ON)
option(WITH_OPENMP          "Enable OpenMP multi-threading"         ON)
option(WITH_DIV32           "Use 32-bit division instead of 64-bit division whenever possible" ON)
option(WITH_MSVC_CRT_STATIC "Link primecount.lib with /MT instead of the default /MD" OFF)
option(WITH_FLOAT128        "Use __float128 (requires libquadmath), increases precision of Li(x) & RiemannR" OFF)
option(WITH_JEMALLOC        "Use jemalloc allocator"                OFF)
```

# Packaging primecount

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

# Man page regeneration

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
