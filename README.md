# primecount

[![Build Status](https://ci.appveyor.com/api/projects/status/github/kimwalisch/primecount?branch=master&svg=true)](https://ci.appveyor.com/project/kimwalisch/primecount)
[![Github Releases](https://img.shields.io/github/release/kimwalisch/primecount.svg)](https://github.com/kimwalisch/primecount/releases)

primecount is a command-line program and [C/C++ library](doc/libprimecount.md)
that counts the primes below an integer x&nbsp;≤&nbsp;10<sup>31</sup> using
**highly optimized** implementations of the combinatorial
[prime counting algorithms](https://en.wikipedia.org/wiki/Prime-counting_function#Algorithms_for_evaluating_%CF%80(x)).

primecount includes implementations of all important combinatorial
prime counting algorithms known up to this date all of which have
been parallelized using [OpenMP](https://en.wikipedia.org/wiki/OpenMP).
primecount contains the first ever open source implementations of
the Deleglise-Rivat algorithm and Xavier Gourdon's algorithm (that works).
primecount also features a [novel load balancer](https://github.com/kimwalisch/primecount/blob/master/src/LoadBalancerS2.cpp)
that is shared amongst all implementations and that scales up to
hundreds of CPU cores. primecount has already been used to compute
several prime counting function [world records](doc/Records.md).

## Installation

The primecount command-line program is available in a few package managers.
For doing development with libprimecount you may need to install
```libprimecount-dev``` or ```libprimecount-devel```.

<table>
    <tr>
        <td><b>Windows:</b></td>
        <td><code>winget install primecount</code></td>
    </tr>
    <tr>
        <td><b>macOS:</b></td>
        <td><code>brew tap kimwalisch/primecount</code><br/>
            <code>brew install primecount</code></td>
    </tr>
    <tr>
        <td><b>Arch Linux:</b></td>
        <td><code>sudo pacman -S primecount</code></td>
    </tr>
    <tr>
        <td><b>Fedora:</b></td>
        <td><code>sudo dnf install primecount</code></td>
    </tr>
    <tr>
        <td><b>FreeBSD:</b></td>
        <td><code>pkg install primecount</code></td>
    </tr>
    <tr>
        <td><b>openSUSE:</b></td>
        <td><code>sudo zypper install primecount</code></td>
    </tr>
</table>

## Build instructions

You need to have installed a C++ compiler and CMake. Ideally
primecount should be compiled using GCC or Clang as these compilers
support both OpenMP (multi-threading library) and 128-bit integers.

```sh
cmake .
make -j
sudo make install
sudo ldconfig
```

* [Detailed build instructions](doc/BUILD.md)

## Usage examples

```sh
# Count the primes below 10^14
primecount 1e14

# Print progress and status information during computation
primecount 1e20 --status

# Count primes using Meissel's algorithm
primecount 2**32 --meissel

# Find the 10^14th prime using 4 threads
primecount 1e14 --nth-prime --threads=4 --time
```

## Command-line options

```
Usage: primecount x [options]
Count the number of primes less than or equal to x (<= 10^31).

Options:

  -d, --deleglise-rivat  Count primes using the Deleglise-Rivat algorithm
  -g, --gourdon          Count primes using Xavier Gourdon's algorithm.
                         This is the default algorithm.
  -l, --legendre         Count primes using Legendre's formula
      --lehmer           Count primes using Lehmer's formula
      --lmo              Count primes using Lagarias-Miller-Odlyzko
  -m, --meissel          Count primes using Meissel's formula
      --Li               Approximate pi(x) using the logarithmic integral
      --Li-inverse       Approximate the nth prime using Li^-1(x)
  -n, --nth-prime        Calculate the nth prime
  -p, --primesieve       Count primes using the sieve of Eratosthenes
      --phi <X> <A>      phi(x, a) counts the numbers <= x that are not
                         divisible by any of the first a primes
      --Ri               Approximate pi(x) using Riemann R
      --Ri-inverse       Approximate the nth prime using Ri^-1(x)
  -s, --status[=NUM]     Show computation progress 1%, 2%, 3%, ...
                         Set digits after decimal point: -s1 prints 99.9%
      --test             Run various correctness tests and exit
      --time             Print the time elapsed in seconds
  -t, --threads=NUM      Set the number of threads, 1 <= NUM <= CPU cores.
                         By default primecount uses all available CPU cores.
  -v, --version          Print version and license information
  -h, --help             Print this help menu
```

<details>
<summary>Advanced options</summary>

```
Advanced options for the Deleglise-Rivat algorithm:

  -a, --alpha=NUM        Set tuning factor: y = x^(1/3) * alpha
      --P2               Compute the 2nd partial sieve function
      --S1               Compute the ordinary leaves
      --S2-trivial       Compute the trivial special leaves
      --S2-easy          Compute the easy special leaves
      --S2-hard          Compute the hard special leaves

Advanced options for Xavier Gourdon's algorithm:

      --alpha-y=NUM      Set tuning factor: y = x^(1/3) * alpha_y
      --alpha-z=NUM      Set tuning factor: z = y * alpha_z
      --AC               Compute the A + C formulas
      --B                Compute the B formula
      --D                Compute the D formula
      --Phi0             Compute the Phi0 formula
      --Sigma            Compute the 7 Sigma formulas
```

</details>

## Benchmarks

<table>
  <tr align="center">
    <td><b>x</b></td>
    <td><b>Prime Count</b></td>
    <td><b>Legendre</b></td>
    <td><b>Meissel</b></td>
    <td><b>Lagarias<br/>Miller<br/>Odlyzko</b></td>
    <td><b>Deleglise<br/>Rivat</b></td>
    <td><b>Gourdon</b></td>
  </tr>
  <tr align="right">
    <td>10<sup>10</sup></td>
    <td>455,052,511</td>
    <td>0.01s</td>
    <td>0.01s</td>
    <td>0.01s</td>
    <td>0.01s</td>
    <td>0.00s</td>
  </tr>
  <tr align="right">
    <td>10<sup>11</sup></td>
    <td>4,118,054,813</td>
    <td>0.01s</td>
    <td>0.01s</td>
    <td>0.01s</td>
    <td>0.01s</td>
    <td>0.01s</td>
  </tr>
  <tr align="right">
    <td>10<sup>12</sup></td>
    <td>37,607,912,018</td>
    <td>0.03s</td>
    <td>0.02s</td>
    <td>0.02s</td>
    <td>0.01s</td>
    <td>0.01s</td>
  </tr>
  <tr align="right">
    <td>10<sup>13</sup></td>
    <td>346,065,536,839</td>
    <td>0.09s</td>
    <td>0.06s</td>
    <td>0.03s</td>
    <td>0.02s</td>
    <td>0.03s</td>
  </tr>
  <tr align="right">
    <td>10<sup>14</sup></td>
    <td>3,204,941,750,802</td>
    <td>0.44s</td>
    <td>0.20s</td>
    <td>0.08s</td>
    <td>0.08s</td>
    <td>0.04s</td>
  </tr>
  <tr align="right">
    <td>10<sup>15</sup></td>
    <td>29,844,570,422,669</td>
    <td>2.33s</td>
    <td>0.89s</td>
    <td>0.29s</td>
    <td>0.16s</td>
    <td>0.11s</td>
  </tr>
  <tr align="right">
    <td>10<sup>16</sup></td>
    <td>279,238,341,033,925</td>
    <td>15.49s</td>
    <td>5.10s</td>
    <td>1.26s</td>
    <td>0.58s</td>
    <td>0.38s</td>
  </tr>
  <tr align="right">
    <td>10<sup>17</sup></td>
    <td>2,623,557,157,654,233</td>
    <td>127.10s</td>
    <td>39.39s</td>
    <td>5.62s</td>
    <td>2.26s</td>
    <td>1.34s</td>
  </tr>
  <tr align="right">
    <td>10<sup>18</sup></td>
    <td>24,739,954,287,740,860</td>
    <td>1,071.14s</td>
    <td>366.93s</td>
    <td>27.19s</td>
    <td>9.96s</td>
    <td>5.35s</td>
  </tr>
  <tr align="right">
    <td>10<sup>19</sup></td>
    <td>234,057,667,276,344,607</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>40.93s</td>
    <td>20.16s</td>
  </tr>
  <tr align="right">
    <td>10<sup>20</sup></td>
    <td>2,220,819,602,560,918,840</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>167.64s</td>
    <td>81.98s</td>
  </tr>
  <tr align="right">
    <td>10<sup>21</sup></td>
    <td>21,127,269,486,018,731,928</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>706.70s</td>
    <td>353.01s</td>
  </tr>
  <tr align="right">
    <td>10<sup>22</sup></td>
    <td>201,467,286,689,315,906,290</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>3,012.10s</td>
    <td>1,350.47s</td>
  </tr>
</table>

The benchmarks above were run on an AMD 7R32 CPU (from 2020) with 16 cores/32 threads
clocked at 3.30GHz. Note that Jan Büthe mentions in <a href="doc/References.md">[11]</a>
that he computed pi(10<sup>25</sup>) in 40,000 CPU core hours using the analytic prime
counting function algorithm. Büthe also mentions that by using additional zeros of the
zeta function the runtime could have potentially been reduced to 4,000 CPU core hours.
However using primecount and Xavier Gourdon's algorithm pi(10<sup>25</sup>) can be computed
in only 460 CPU core hours on an AMD Ryzen 3950X CPU!

## Performance tips

If you have an x64 CPU and you have installed primecount using the package manager of
your Linux distribution, then it is possible that the ```POPCNT``` instruction has been
disabled in order to ensure that primecount works on very old CPUs. Unfortunately this
decreases performance by up to 50%. On the other hand, if you compile primecount from
source the ```POPCNT``` instruction will be enabled by default. The fastest primecount
binary can be built using the Clang compiler and the ```-march=native``` option.

```bash
CXX=clang++ CXXFLAGS="-march=native" cmake .
make -j
```

By default primecount scales nicely up until 10<sup>20</sup> on current x64 CPUs.
For larger values primecount's large memory usage causes many
[TLB (translation lookaside buffer)](https://en.wikipedia.org/wiki/Translation_lookaside_buffer)
cache misses that significantly deteriorate primecount's performance.
Fortunately the Linux kernel allows to enable
[transparent huge pages](https://www.kernel.org/doc/html/latest/admin-guide/mm/transhuge.html)
so that large memory allocations will automatically be done using huge
pages instead of ordinary pages which dramatically reduces the number of
TLB cache misses.

```bash
# Enable transparent huge pages until next reboot
sudo bash -c 'echo always > /sys/kernel/mm/transparent_hugepage/enabled'
```

## Algorithms

<table>
  <tr>
    <td>Legendre's Formula</td>
    <td><img src="https://kimwalisch.github.io/primecount/formulas/pi_legendre.svg" height="20" align="absmiddle"/></td>
  </tr>
  <tr>
    <td>Meissel's Formula</td>
    <td><img src="https://kimwalisch.github.io/primecount/formulas/pi_meissel.svg" height="20" align="absmiddle"/></td>
  </tr>
  <tr>
    <td>Lehmer's Formula</td>
    <td><img src="https://kimwalisch.github.io/primecount/formulas/pi_lehmer.svg" height="20" align="absmiddle"/></td>
  </tr>
  <tr>
    <td>LMO Formula</td>
    <td><img src="https://kimwalisch.github.io/primecount/formulas/pi_lmo.svg" height="20" align="absmiddle"/></td>
  </tr>
</table>

Up until the early 19th century the most efficient known method for
counting primes was the sieve of Eratosthenes which has a running time of
<img src="https://kimwalisch.github.io/primecount/formulas/Oxloglogx.svg" height="20" align="absmiddle"/>
operations. The first improvement to this bound was Legendre's formula
(1830) which uses the inclusion-exclusion principle to calculate the
number of primes below x without enumerating the individual primes.
Legendre's formula has a running time of
<img src="https://kimwalisch.github.io/primecount/formulas/Ox.svg" height="20" align="absmiddle"/>
operations and uses
<img src="http://kimwalisch.github.io/primecount/formulas/Osqrtx.svg" height="20" align="absmiddle"/>
space. In 1870 E. D. F. Meissel improved Legendre's formula by setting
<img src="https://kimwalisch.github.io/primecount/formulas/apisqrt3x.svg" height="20" align="absmiddle"/>
and by adding the correction term
<img src="https://kimwalisch.github.io/primecount/formulas/P2xa.svg" height="20" align="absmiddle"/>.
Meissel's formula has a running time of
<img src="https://kimwalisch.github.io/primecount/formulas/Omeissel.svg" height="20" align="absmiddle"/>
operations and uses
<img src="https://kimwalisch.github.io/primecount/formulas/Osqrtxlogx.svg" height="20" align="absmiddle"/>
space. In 1959 D. H. Lehmer extended Meissel's formula and slightly improved the running time to
<img src="https://kimwalisch.github.io/primecount/formulas/Olehmer.svg" height="20" align="absmiddle"/>
operations and
<img src="https://kimwalisch.github.io/primecount/formulas/Osqrtxlogx.svg" height="20" align="absmiddle"/>
space. In 1985 J. C. Lagarias, V. S. Miller and A. M. Odlyzko published a new
algorithm based on Meissel's formula which has a lower runtime complexity of
<img src="https://kimwalisch.github.io/primecount/formulas/Oroot23xlogx.svg" height="20" align="absmiddle"/>
operations and which uses only
<img src="https://kimwalisch.github.io/primecount/formulas/Osqrt3xlog2x.svg" height="20" align="absmiddle"/>
space.

primecount's Legendre, Meissel and Lehmer implementations are based
on Hans Riesel's book <a href="doc/References.md">[5]</a>,
its Lagarias-Miller-Odlyzko and Deleglise-Rivat implementations are
based on Tomás Oliveira's paper <a href="doc/References.md">[9]</a>
and the implementation of Xavier Gourdon's algorithm is based
on Xavier Gourdon's paper <a href="doc/References.md">[7]</a>.
primecount's implementation of the so-called hard special leaves is different
from the algorithms that have been described in any of the combinatorial
prime counting papers so far. Instead of using a binary indexed tree
for counting which is very cache inefficient primecount uses a linear
counter array in combination with the POPCNT instruction which is more
cache efficient and much faster. The
[Hard-Special-Leaves.md](doc/Hard-Special-Leaves.md) document contains more
information. primecount's [easy special leaf](doc/Easy-Special-Leaves.md)
implementation and its [partial sieve function](doc/Partial-Sieve-Function.md)
implementation also contain significant improvements.

## Fast nth prime calculation

The most efficient known method for calculating the nth prime is a
combination of the prime counting function and a prime sieve. The idea
is to closely approximate the nth prime (e.g. using the inverse
logarithmic integral
<img src="https://kimwalisch.github.io/primecount/formulas/Li-1n.svg" height="20" align="absmiddle"/>
or the inverse Riemann R function
<img src="https://kimwalisch.github.io/primecount/formulas/RiemannR-1.svgz" height="20" align="absmiddle"/>)
and then count the primes up to this guess using the prime counting
function. Once this is done one starts sieving (e.g. using the
segmented sieve of Eratosthenes) from there on until one finds the
actual nth prime. The author has implemented ```primecount::nth_prime(n)```
this way (option: ```--nth-prime```), it finds the nth prime in
<img src="https://kimwalisch.github.io/primecount/formulas/Oroot23xlog2x.svg" height="20" align="absmiddle"/>
operations using
<img src="https://kimwalisch.github.io/primecount/formulas/Opisqrtx.svg" height="20" align="absmiddle"/>
space.

## C API

Include the ```<primecount.h>``` header to use primecount's C API.
All functions that are part of primecount's C API return ```-1``` in case an
error occurs and print the corresponding error message to the standard error
stream.

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

* [C API documentation](doc/libprimecount.md#libprimecount)
* [libprimecount build instructions](doc/libprimecount.md#build-instructions)

## C++ API

Include the ```<primecount.hpp>``` header to use primecount's C++ API.
All functions that are part of primecount's C++ API throw a
```primecount_error``` exception (which is derived from
```std::exception```) in case an error occurs.

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

* [C++ API documentation](doc/libprimecount.md#libprimecount)
* [libprimecount build instructions](doc/libprimecount.md#build-instructions)

## Bindings for other languages

primesieve natively supports C and C++ and has bindings available for:

<table>
    <tr>
        <td><b>Common Lisp:</b></td>
        <td><a href="https://github.com/AaronChen0/cl-primecount">cl-primecount</a></td>
    </tr>
    <tr>
        <td><b>Julia:</b></td>
        <td><a href="https://github.com/JuliaBinaryWrappers/primecount_jll.jl">primecount_jll.jl</a></td>
    </tr>
    <tr>
        <td><b>Haskell:</b></td>
        <td><a href="https://github.com/pgujjula/primecount-haskell">primecount-haskell</a></td>
    </tr>
    <tr>
        <td><b>Python:</b></td>
        <td><a href="https://github.com/dimpase/primecountpy">primecountpy</a></td>
    </tr>
    <tr>
        <td><b>Python:</b></td>
        <td><a href="https://github.com/hearot/primecount-python">primecount-python</a></td>
    </tr>
    <tr>
        <td><b>Rust:</b></td>
        <td><a href="https://github.com/maitbayev/primecount-rs">primecount-rs</a></td>
    </tr>
</table>

Many thanks to the developers of these bindings!
