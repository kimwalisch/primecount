# primecount

[![Build Status](https://travis-ci.org/kimwalisch/primecount.svg)](https://travis-ci.org/kimwalisch/primecount)
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
primecount also features a [novel load balancer](https://github.com/kimwalisch/primecount/blob/master/src/LoadBalancer.cpp)
that is shared amongst all implementations and that scales up to
hundreds of CPU cores. primecount has already been used to compute
several [world records](doc/Records.md).

## Build instructions

You need to have installed a C++ compiler and CMake. Ideally
primecount should be compiled using a C++ compiler that supports both
OpenMP and 128-bit integers (e.g. GCC, Clang, Intel C++ Compiler).

```sh
cmake .
make -j
sudo make install
```

* [Detailed build instructions](doc/BUILD.md)

## Binaries

Below are the latest precompiled primecount binaries for
Windows, Linux and macOS. These binaries are statically linked
and require a CPU which supports the POPCNT instruction (2008 or
later).

* [primecount-6.0-win64.zip](https://github.com/kimwalisch/primecount/releases/download/v6.0/primecount-6.0-win64.zip), 589 kB
* [primecount-6.0-linux-x64.tar.xz](https://github.com/kimwalisch/primecount/releases/download/v6.0/primecount-6.0-linux-x64.tar.xz), 800 kB
* [primecount-6.0-macOS-x64.zip](https://github.com/kimwalisch/primecount/releases/download/v6.0/primecount-6.0-macOS-x64.zip), 392 kB
* Binaries with backup functionality are available [here](https://github.com/kimwalisch/primecount/tree/backup3#primecount-backup)

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
    <td>0.02s</td>
    <td>0.01s</td>
    <td>0.00s</td>
    <td>0.00s</td>
    <td>0.00s</td>
  </tr>
  <tr align="right">
    <td>10<sup>11</sup></td>
    <td>4,118,054,813</td>
    <td>0.02s</td>
    <td>0.01s</td>
    <td>0.01s</td>
    <td>0.01s</td>
    <td>0.01s</td>
  </tr>
  <tr align="right">
    <td>10<sup>12</sup></td>
    <td>37,607,912,018</td>
    <td>0.09s</td>
    <td>0.05s</td>
    <td>0.02s</td>
    <td>0.01s</td>
    <td>0.01s</td>
  </tr>
  <tr align="right">
    <td>10<sup>13</sup></td>
    <td>346,065,536,839</td>
    <td>0.39s</td>
    <td>0.19s</td>
    <td>0.05s</td>
    <td>0.03s</td>
    <td>0.02s</td>
  </tr>
  <tr align="right">
    <td>10<sup>14</sup></td>
    <td>3,204,941,750,802</td>
    <td>2.24s</td>
    <td>0.96s</td>
    <td>0.16s</td>
    <td>0.12s</td>
    <td>0.07s</td>
  </tr>
  <tr align="right">
    <td>10<sup>15</sup></td>
    <td>29,844,570,422,669</td>
    <td>14.58s</td>
    <td>5.96s</td>
    <td>0.63s</td>
    <td>0.38s</td>
    <td>0.21s</td>
  </tr>
  <tr align="right">
    <td>10<sup>16</sup></td>
    <td>279,238,341,033,925</td>
    <td>111.81s</td>
    <td>42.91s</td>
    <td>2.83s</td>
    <td>1.59s</td>
    <td>0.76s</td>
  </tr>
  <tr align="right">
    <td>10<sup>17</sup></td>
    <td>2,623,557,157,654,233</td>
    <td>938.07s</td>
    <td>352.39s</td>
    <td>12.70s</td>
    <td>4.72s</td>
    <td>2.83s</td>
  </tr>
  <tr align="right">
    <td>10<sup>18</sup></td>
    <td>24,739,954,287,740,860</td>
    <td>8,268.52s</td>
    <td>3,144.72s</td>
    <td>57.31s</td>
    <td>18.96s</td>
    <td>10.77s</td>
  </tr>
  <tr align="right">
    <td>10<sup>19</sup></td>
    <td>234,057,667,276,344,607</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>89.02s</td>
    <td>45.56s</td>
  </tr>
  <tr align="right">
    <td>10<sup>20</sup></td>
    <td>2,220,819,602,560,918,840</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>385.61s</td>
    <td>188.70s</td>
  </tr>
  <tr align="right">
    <td>10<sup>21</sup></td>
    <td>21,127,269,486,018,731,928</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>1,652.66s</td>
    <td>800.57s</td>
  </tr>
  <tr align="right">
    <td>10<sup>22</sup></td>
    <td>201,467,286,689,315,906,290</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>7,355.69s</td>
    <td>3,260.95s</td>
  </tr>
</table>

The benchmarks above were run on a system with an Intel Xeon Platinum 8275CL CPU
from 2019 using 8 CPU cores (16 threads) clocked at 3.00GHz. Note that Jan Büthe
mentions in <a href="doc/References.md">[11]</a> that he computed pi(10<sup>25</sup>)
in 40,000 CPU core hours using the analytic prime counting function algorithm.
Büthe also mentions that by using additional zeros of the zeta function the runtime
could have potentially been reduced to 4,000 CPU core hours. However using primecount
and Xavier Gourdon's algorithm pi(10<sup>25</sup>) can be computed in only 460 CPU
core hours on an AMD Ryzen 3950X CPU!

## Performance tips

primecount uses the OpenMP multi-threading library. In OpenMP waiting threads are
usually busy-waiting for a short amount of time (by spinning) before being put to sleep.
This setting can be altered using the ```OMP_WAIT_POLICY``` environment variable. For
primecount it is best to set ```OMP_WAIT_POLICY``` to ```PASSIVE``` in order to prevent
the threads from busy waiting. This setting can provide a large speedup for small
to medium sized computations below 10<sup>22</sup>.

```bash
export OMP_WAIT_POLICY=PASSIVE
```

By default primecount scales nicely up until 10<sup>23</sup> on current x64 CPUs.
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
primecount's implementation of the special leaves formula is different
from the algorithms that have been described in any of the combinatorial
prime counting papers so far. Instead of using a binary indexed tree
for counting which is very cache inefficient primecount uses a linear
counters array in combination with the POPCNT instruction which is more
cache efficient and much faster. The
[Special-Leaves.md](doc/Special-Leaves.md) document contains more
information.

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

* [C API reference](doc/libprimecount.md#c-example)
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

* [C++ API reference](doc/libprimecount.md#c-example-1)
* [libprimecount build instructions](doc/libprimecount.md#build-instructions)
