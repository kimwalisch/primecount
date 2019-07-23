# primecount

[![Build Status](https://travis-ci.org/kimwalisch/primecount.svg)](https://travis-ci.org/kimwalisch/primecount)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/kimwalisch/primecount?branch=master&svg=true)](https://ci.appveyor.com/project/kimwalisch/primecount)
[![Github Releases](https://img.shields.io/github/release/kimwalisch/primecount.svg)](https://github.com/kimwalisch/primecount/releases)

primecount is a command-line program and [C++ library](doc/libprimecount.md)
that counts the primes below an integer x&nbsp;≤&nbsp;10<sup>31</sup> using
**highly optimized** implementations of the
[prime counting function](https://en.wikipedia.org/wiki/Prime-counting_function)
(combinatorial methods). primecount includes implementations of the
algorithms of Legendre, Meissel, Lehmer, Lagarias-Miller-Odlyzko, Deleglise-Rivat
and Gourdon all of which have been parallelized using
[OpenMP](https://en.wikipedia.org/wiki/OpenMP). The Deleglise-Rivat
implementation has also been 
[distributed](https://github.com/kimwalisch/primecount/blob/master/doc/primecount-MPI.md#primecount-mpi)
using MPI.

primecount contains the **first ever** parallel open source
implementation of the Deleglise-Rivat algorithm and it features a
[novel load balancer](https://github.com/kimwalisch/primecount/blob/master/src/LoadBalancer.cpp)
which scales up to hundreds of CPU cores. primecount has already been
used to compute several world records e.g.
[pi(10<sup>27</sup>)](http://www.mersenneforum.org/showthread.php?t=20473) and
[nth_prime(10<sup>24</sup>)](https://oeis.org/A006988), more will
hopefully follow!

## Build instructions

You need to have installed a C++ compiler and CMake. Ideally
primecount should be compiled using a C++ compiler that supports both
OpenMP and 128-bit integers (e.g. GCC, Clang, Intel C++ Compiler).

```sh
cmake .
make -j
sudo make install
```

To build primecount using
[MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface)
support for distributing computations onto cluster nodes use:

```sh
cmake -DWITH_MPI=ON .
```

[primecount-MPI.md](doc/primecount-MPI.md) contains more information.

## Binaries

Below are the latest precompiled primecount binaries for
Windows, Linux and macOS. These binaries are statically linked
and require a CPU which supports the POPCNT instruction (2010 or
later).

* [primecount-4.8-win64.zip](https://github.com/kimwalisch/primecount/releases/download/v4.8/primecount-4.8-win64.zip), 510 KB
* [primecount-4.8-linux-x64.tar.xz](https://github.com/kimwalisch/primecount/releases/download/v4.8/primecount-4.8-linux-x64.tar.xz), 824 KB
* [primecount-4.8-macOS-x64.zip](https://github.com/kimwalisch/primecount/releases/download/v4.8/primecount-4.8-macOS-x64.zip), 346 KB
* Binaries with backup functionality are available [here](https://github.com/kimwalisch/primecount/tree/backup2#primecount-backup)

## Usage examples

Open a terminal and run primecount using e.g.:
```sh
# Count the primes below 10^14
./primecount 1e14

# Print progress and status information during computation
./primecount 1e20 --status

# Count primes using Meissel's algorithm
./primecount 2**32 --meissel

# Find the 10^14th prime using 4 threads
./primecount 1e14 --nthprime --threads=4 --time
```

## Command-line options

```
Usage: primecount x [OPTION]...
Count the primes below x <= 10^31 using fast implementations of the
combinatorial prime counting function.

Options:

  -d,    --deleglise_rivat  Count primes using Deleglise-Rivat algorithm
  -g,    --gourdon          Count primes using Xavier Gourdon's algorithm
         --legendre         Count primes using Legendre's formula
         --lehmer           Count primes using Lehmer's formula
  -l,    --lmo              Count primes using Lagarias-Miller-Odlyzko
  -m,    --meissel          Count primes using Meissel's formula
         --Li               Approximate pi(x) using the logarithmic integral
         --Li_inverse       Approximate the nth prime using Li^-1(x)
  -n,    --nthprime         Calculate the nth prime
  -p,    --primesieve       Count primes using the sieve of Eratosthenes
         --phi=<a>          phi(x, a) counts the numbers <= x that are
                            not divisible by any of the first a primes
         --Ri               Approximate pi(x) using Riemann R
         --Ri_inverse       Approximate the nth prime using Ri^-1(x)
  -s[N], --status[=N]       Show computation progress 1%, 2%, 3%, ...
                            [N] digits after decimal point e.g. N=1, 99.9%
         --test             Run various correctness tests and exit
         --time             Print the time elapsed in seconds
  -t<N>, --threads=<N>      Set the number of threads, 1 <= N <= CPU cores
  -v,    --version          Print version and license information
  -h,    --help             Print this help menu

Advanced options for the Deleglise-Rivat algorithm:

  -a<N>, --alpha=<N>        Tuning factor, 1 <= alpha <= x^(1/6)
         --P2               Only compute the 2nd partial sieve function
         --S1               Only compute the ordinary leaves
         --S2_trivial       Only compute the trivial special leaves
         --S2_easy          Only compute the easy special leaves
         --S2_hard          Only compute the hard special leaves

Advanced options for Xavier Gourdon's algorithm:

         --alpha_y=<N>      Tuning factor, 1 <= alpha_y <= x^(1/6)
                            with y = x^(1/3) * alpha_y
         --alpha_z=<N>      Tuning factor, 1 <= alpha_z <= alpha_y
                            with z = y * alpha_z
```

## Benchmarks

<table>
  <tr align="center">
    <td><b>x</b></td>
    <td><b>Prime Count</b></td>
    <td><b>Legendre</b></td>
    <td><b>Meissel</b></td>
    <td><b>Lagarias<br/>Miller<br/>Odlyzko</b></td>
    <td><b>Deleglise<br/>Rivat</b></td>
  </tr>
  <tr align="right">
    <td>10<sup>10</sup></td>
    <td>455,052,511</td>
    <td>0.02s</td>
    <td>0.01s</td>
    <td>0.01s</td>
    <td>0.01s</td>
  </tr>
  <tr align="right">
    <td>10<sup>11</sup></td>
    <td>4,118,054,813</td>
    <td>0.03s</td>
    <td>0.02s</td>
    <td>0.01s</td>
    <td>0.01s</td>
  </tr>
  <tr align="right">
    <td>10<sup>12</sup></td>
    <td>37,607,912,018</td>
    <td>0.10s</td>
    <td>0.06s</td>
    <td>0.03s</td>
    <td>0.02s</td>
  </tr>
  <tr align="right">
    <td>10<sup>13</sup></td>
    <td>346,065,536,839</td>
    <td>0.46s</td>
    <td>0.26s</td>
    <td>0.06s</td>
    <td>0.04s</td>
  </tr>
  <tr align="right">
    <td>10<sup>14</sup></td>
    <td>3,204,941,750,802</td>
    <td>2.66s</td>
    <td>1.28s</td>
    <td>0.19s</td>
    <td>0.10s</td>
  </tr>
  <tr align="right">
    <td>10<sup>15</sup></td>
    <td>29,844,570,422,669</td>
    <td>17.51s</td>
    <td>7.65s</td>
    <td>0.73s</td>
    <td>0.35s</td>
  </tr>
  <tr align="right">
    <td>10<sup>16</sup></td>
    <td>279,238,341,033,925</td>
    <td>136.16s</td>
    <td>55.30s</td>
    <td>2.95s</td>
    <td>1.26s</td>
  </tr>
  <tr align="right">
    <td>10<sup>17</sup></td>
    <td>2,623,557,157,654,233</td>
    <td>1,155.34s</td>
    <td>456.08s</td>
    <td>12.49s</td>
    <td>4.57s</td>
  </tr>
  <tr align="right">
    <td>10<sup>18</sup></td>
    <td>24,739,954,287,740,860</td>
    <td>10,239.22s</td>
    <td>3,974.21s</td>
    <td>55.87s</td>
    <td>18.97s</td>
  </tr>
  <tr align="right">
    <td>10<sup>19</sup></td>
    <td>234,057,667,276,344,607</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>86.13s</td>
  </tr>
  <tr align="right">
    <td>10<sup>20</sup></td>
    <td>2,220,819,602,560,918,840</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>373.09s</td>
  </tr>
  <tr align="right">
    <td>10<sup>21</sup></td>
    <td>21,127,269,486,018,731,928</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>1,747.57s</td>
  </tr>
  <tr align="right">
    <td>10<sup>22</sup></td>
    <td>201,467,286,689,315,906,290</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>8,215.59s</td>
  </tr>
</table>

The benchmarks above were run on an AMD Ryzen 7 1700 CPU (8 x 3.0
GHz, 3.7 GHz Turbo) from 2017 and primecount was compiled using
GCC 6.3. Only the Deleglise-Rivat implementation supports
numbers&nbsp;>&nbsp;2<sup>63</sup>.

## Performance tips

primecount scales nicely up until 10<sup>24</sup> on current CPUs. For larger
values primecount's large memory usage causes many
[TLB (translation lookaside buffer)](https://en.wikipedia.org/wiki/Translation_lookaside_buffer)
cache misses that severely deteriorate primecount's performance.
Fortunately the Linux kernel allows to enable
[transparent huge pages](https://www.kernel.org/doc/html/latest/admin-guide/mm/transhuge.html)
so that large memory allocations will automatically be done using huge
pages instead of ordinary pages which dramatically reduces the number of
TLB cache misses. I have measured a speed up of more than 2x for
S1(10^26) when transparent huge pages are enabled.

```bash
sudo su

# Enable transparent huge pages until next reboot
echo always > /sys/kernel/mm/transparent_hugepage/enabled
echo always > /sys/kernel/mm/transparent_hugepage/defrag
```

## C++ library

primecount can be built as a static and shared C++ library for use in
other math projects.

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

[libprimecount.md](doc/libprimecount.md) contains more information.

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

<p>Up until the early 19th century the most efficient known method for
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
space.</p>
<p>primecount's Legendre, Meissel and Lehmer implementations are based
on Hans Riesel's book <a href="doc/References.md">[5]</a>,
its Lagarias-Miller-Odlyzko and Deleglise-Rivat implementations are
based on Tomás Oliveira's paper <a href="doc/References.md">[9]</a>.</p>

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
this way (option: ```--nthprime```), it finds the nth prime in
<img src="https://kimwalisch.github.io/primecount/formulas/Oroot23xlog2x.svg" height="20" align="absmiddle"/>
operations using
<img src="https://kimwalisch.github.io/primecount/formulas/Opisqrtx.svg" height="20" align="absmiddle"/>
space.
