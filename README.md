# primecount

[![Build status](https://github.com/kimwalisch/primecount/actions/workflows/ci.yml/badge.svg)](https://github.com/kimwalisch/primecount/actions/workflows/ci.yml) [![Build status](https://github.com/kimwalisch/primecount/actions/workflows/benchmark.yml/badge.svg)](https://github.com/kimwalisch/primecount/actions/workflows/benchmark.yml)
[![Github Releases](https://img.shields.io/github/release/kimwalisch/primecount.svg)](https://github.com/kimwalisch/primecount/releases)
[![C API Documentation](https://img.shields.io/badge/docs-C_API-blue)](doc/libprimecount.md)
[![C++ API Documentation](https://img.shields.io/badge/docs-C++_API-blue)](doc/libprimecount.md)

primecount is a command-line program and C/C++ library that counts the number of
primes&nbsp;≤&nbsp;x (maximum 10<sup>31</sup>) using **highly optimized** implementations of the
[combinatorial prime counting algorithms](https://en.wikipedia.org/wiki/Prime-counting_function#Algorithms_for_evaluating_%CF%80(x)).

primecount includes implementations of all important combinatorial prime counting algorithms
known up to this date all of which have been parallelized using
[OpenMP](https://en.wikipedia.org/wiki/OpenMP). primecount contains the first ever open
source implementations of the Deleglise-Rivat algorithm and Xavier Gourdon's algorithm (that works).
primecount also features a [novel load balancer](https://github.com/kimwalisch/primecount/blob/master/src/LoadBalancerS2.cpp)
that is shared amongst all implementations and that scales up to hundreds of CPU cores. primecount
has already been used to compute several prime counting function [world records](doc/Records.md)!

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
        <td><code>brew install primecount</code></td>
    </tr>
    <tr>
        <td><b>Arch Linux:</b></td>
        <td><code>sudo pacman -S primecount</code></td>
    </tr>
    <tr>
        <td><b>Debian/Ubuntu:</b></td>
        <td><code>sudo apt install primecount</code></td>
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
cmake --build . --parallel
sudo cmake --install .
sudo ldconfig
```

* [Detailed build instructions](doc/BUILD.md)

## Usage examples

```sh
# Count the primes ≤ 10^14
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

  -d, --deleglise-rivat    Count primes using the Deleglise-Rivat algorithm
      --double-check       Recompute pi(x) with alternative alpha tuning
                           factor(s) to verify the first result.
  -g, --gourdon            Count primes using Xavier Gourdon's algorithm.
                           This is the default algorithm.
  -l, --legendre           Count primes using Legendre's formula
      --lehmer             Count primes using Lehmer's formula
      --lmo                Count primes using Lagarias-Miller-Odlyzko
  -m, --meissel            Count primes using Meissel's formula
      --Li                 Eulerian logarithmic integral function
      --Li-inverse         Approximate the nth prime using Li^-1(x)
  -n, --nth-prime          Calculate the nth prime
  -p, --primesieve         Count primes using the sieve of Eratosthenes
      --phi <X> <A>        phi(x, a) counts the numbers <= x that are not
                           divisible by any of the first a primes
  -R, --RiemannR           Approximate pi(x) using the Riemann R function
      --RiemannR-inverse   Approximate the nth prime using R^-1(x)
  -s, --status[=NUM]       Show computation progress 1%, 2%, 3%, ...
                           Set digits after decimal point: -s1 prints 99.9%
      --test               Run various correctness tests and exit
      --time               Print the time elapsed in seconds
  -t, --threads=NUM        Set the number of threads, 1 <= NUM <= CPU cores.
                           By default primecount uses all available CPU cores.
  -v, --version            Print version and license information
  -h, --help               Print this help menu
```

<details>
<summary>Advanced options</summary>

```
Advanced options for the Deleglise-Rivat algorithm:

  -a, --alpha=NUM          Set tuning factor: y = x^(1/3) * alpha
      --P2                 Compute the 2nd partial sieve function
      --S1                 Compute the ordinary leaves
      --S2-trivial         Compute the trivial special leaves
      --S2-easy            Compute the easy special leaves
      --S2-hard            Compute the hard special leaves

Advanced options for Xavier Gourdon's algorithm:

      --alpha-y=NUM        Set tuning factor: y = x^(1/3) * alpha_y
      --alpha-z=NUM        Set tuning factor: z = y * alpha_z
      --AC                 Compute the A + C formulas
      --B                  Compute the B formula
      --D                  Compute the D formula
      --Phi0               Compute the Phi0 formula
      --Sigma              Compute the 7 Sigma formulas
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
    <td>0.00s</td>
    <td>0.00s</td>
    <td>0.00s</td>
    <td>0.00s</td>
  </tr>
  <tr align="right">
    <td>10<sup>11</sup></td>
    <td>4,118,054,813</td>
    <td>0.01s</td>
    <td>0.01s</td>
    <td>0.01s</td>
    <td>0.01s</td>
    <td>0.00s</td>
  </tr>
  <tr align="right">
    <td>10<sup>12</sup></td>
    <td>37,607,912,018</td>
    <td>0.02s</td>
    <td>0.01s</td>
    <td>0.01s</td>
    <td>0.01s</td>
    <td>0.01s</td>
  </tr>
  <tr align="right">
    <td>10<sup>13</sup></td>
    <td>346,065,536,839</td>
    <td>0.03s</td>
    <td>0.02s</td>
    <td>0.02s</td>
    <td>0.02s</td>
    <td>0.01s</td>
  </tr>
  <tr align="right">
    <td>10<sup>14</sup></td>
    <td>3,204,941,750,802</td>
    <td>0.11s</td>
    <td>0.05s</td>
    <td>0.03s</td>
    <td>0.03s</td>
    <td>0.02s</td>
  </tr>
  <tr align="right">
    <td>10<sup>15</sup></td>
    <td>29,844,570,422,669</td>
    <td>0.45s</td>
    <td>0.21s</td>
    <td>0.14s</td>
    <td>0.13s</td>
    <td>0.06s</td>
  </tr>
  <tr align="right">
    <td>10<sup>16</sup></td>
    <td>279,238,341,033,925</td>
    <td>3.09s</td>
    <td>1.12s</td>
    <td>0.41s</td>
    <td>0.31s</td>
    <td>0.20s</td>
  </tr>
  <tr align="right">
    <td>10<sup>17</sup></td>
    <td>2,623,557,157,654,233</td>
    <td>25.28s</td>
    <td>8.84s</td>
    <td>1.81s</td>
    <td>1.27s</td>
    <td>0.51s</td>
  </tr>
  <tr align="right">
    <td>10<sup>18</sup></td>
    <td>24,739,954,287,740,860</td>
    <td>214.63s</td>
    <td>78.00s</td>
    <td>8.18s</td>
    <td>5.33s</td>
    <td>2.00s</td>
  </tr>
  <tr align="right">
    <td>10<sup>19</sup></td>
    <td>234,057,667,276,344,607</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>24.40s</td>
    <td>8.12s</td>
  </tr>
  <tr align="right">
    <td>10<sup>20</sup></td>
    <td>2,220,819,602,560,918,840</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>113.60s</td>
    <td>32.87s</td>
  </tr>
  <tr align="right">
    <td>10<sup>21</sup></td>
    <td>21,127,269,486,018,731,928</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>500.51s</td>
    <td>134.21s</td>
  </tr>
  <tr align="right">
    <td>10<sup>22</sup></td>
    <td>201,467,286,689,315,906,290</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>2,198.92s</td>
    <td>552.17s</td>
  </tr>
</table>

The benchmarks above were run on an AMD EPYC Zen4 CPU from 2023 with 32 CPU cores (no Hyper-Threading)
clocked at 3.7 GHz. Note that Jan Büthe mentions in <a href="doc/References.md">[11]</a>
that he computed $\pi(10^{25})$ in 40,000 CPU core hours using the analytic prime
counting function algorithm. Büthe also mentions that by using additional zeros of the
zeta function the runtime could have potentially been reduced to 4,000 CPU core hours.
However using primecount and Xavier Gourdon's algorithm $\pi(10^{25})$ can be computed
in only 357 CPU core hours on the AMD EPYC Zen4 CPU!

## Algorithms

<table>
  <tr>
    <td>Legendre's Formula</td>
    <td>$\pi(x)=\pi(\sqrt{x})+\phi(x,\pi(\sqrt{x}))-1$</td>
  </tr>
  <tr>
    <td>Meissel's Formula</td>
    <td>$\pi(x)=\pi(\sqrt[3]{x})+\phi(x,\pi(\sqrt[3]{x}))-\mathrm{P_2}(x,\pi(\sqrt[3]{x}))-1$</td>
  </tr>
  <tr>
    <td>Lehmer's Formula</td>
    <td>$\pi(x)=\pi(\sqrt[4]{x})+\phi(x,\pi(\sqrt[4]{x}))-\mathrm{P_2}(x,\pi(\sqrt[4]{x}))-\mathrm{P_3}(x,\pi(\sqrt[4]{x}))-1$</td>
  </tr>
  <tr>
    <td>LMO Formula</td>
    <td>$\pi(x)=\pi(\sqrt[3]{x})+\mathrm{S_1}(x,\pi(\sqrt[3]{x}))+\mathrm{S_2}(x,\pi(\sqrt[3]{x}))-\mathrm{P_2}(x,\pi(\sqrt[3]{x}))-1$</td>
  </tr>
</table>

Up until the early 19th century the most efficient known method for counting primes was the
sieve of Eratosthenes which has a running time of $O(x\ \log\ \log\ x)$ operations. The first
improvement to this bound was Legendre's formula (1830) which uses the inclusion-exclusion
principle to calculate the number of primes below x without enumerating the individual
primes. Legendre's formula has a running time of $O(x)$ operations and uses $O(\sqrt{x}/\log{x})$
space. In 1870 E. D. F. Meissel improved Legendre's formula by setting $a=\pi(\sqrt[3]{x})$
and by adding the correction term $\mathrm{P_2}(x,a)$, Meissel's formula has a running time
of $O(x/\log^3{x})$ operations and uses $O(\sqrt[3]{x})$ space. In 1959 D. H. Lehmer
extended Meissel's formula and slightly improved the running time to $O(x/\log^4{x})$
operations and $O(x^{\frac{3}{8}})$ space. In 1985 J. C. Lagarias, V. S. Miller and A. M.
Odlyzko published a new algorithm based on Meissel's formula which has a lower runtime
complexity of $O(x^{\frac{2}{3}}/\log{x})$ operations and which uses only
$O(\sqrt[3]{x}\ \log^2{x})$ space.

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
cache efficient and much faster. [Hard-Special-Leaves.pdf](doc/Hard-Special-Leaves.pdf)
explains this algorithm in more detail. The novelties of primecount's easy special leaf
implementation are described in [Easy-Special-Leaves.pdf](doc/Easy-Special-Leaves.pdf),
and its partial sieve function implementation is described in
[Partial-Sieve-Function.pdf](doc/Partial-Sieve-Function.pdf).

## Fast nth prime calculation

The most efficient known method for calculating the nth prime is a combination
of the prime counting function and a prime sieve. The idea is to closely
approximate the nth prime e.g. using the inverse logarithmic integral
$\mathrm{Li}^{-1}(n)$ or the inverse Riemann R function $\mathrm{R}^{-1}(n)$
and then count the primes up to this guess using the prime counting function.
Once this is done one starts sieving (e.g. using the segmented sieve of
Eratosthenes) from there on until one finds the actual nth prime. The author
has implemented ```primecount::nth_prime(n)``` this way
(option: ```--nth-prime```), it finds the nth prime in $O(x^{\frac{2}{3}}/\log^2{x})$
operations using $O(\sqrt[3]{x}\ \log^3{x})$ space.

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
    printf("primes <= 1000: %ld\n", pix);

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
    std::cout << "primes <= 1000: " << pix << std::endl;

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
        <td><b>Lua:</b></td>
        <td><a href="https://github.com/ishandutta2007/lua-primecount">lua-primecount</a></td>
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
        <td><b>Rust:</b></td>
        <td><a href="https://github.com/maitbayev/primecount-rs">primecount-rs</a></td>
    </tr>
</table>

Many thanks to the developers of these bindings!

## Sponsors

Thanks to all current and past [sponsors of primecount](https://github.com/sponsors/kimwalisch)! Your donations help me purchase (or rent) the latest CPUs and ensure primecount runs at maximum performance on them. Your donations also motivate me to continue maintaining primecount.

<a href="https://github.com/AndrewVSutherland"><img src="https://images.weserv.nl/?url=avatars.githubusercontent.com/u/11425002?h=60&w=60&fit=cover&mask=circle"></img></a>
<a href="https://github.com/wolframresearch"><img src="https://images.weserv.nl/?url=avatars.githubusercontent.com/u/11549616?h=60&w=60&fit=cover&mask=circle"></img></a>
<a href="https://github.com/AlgoWin"><img src="https://images.weserv.nl/?url=avatars.githubusercontent.com/u/44401099?h=60&w=60&fit=cover&mask=circle"></img></a>
<a href="https://github.com/sethtroisi"><img src="https://images.weserv.nl/?url=avatars.githubusercontent.com/u/10172976?h=60&w=60&fit=cover&mask=circle"></img></a>
<a href="https://github.com/entersoftone"><img src="https://images.weserv.nl/?url=avatars.githubusercontent.com/u/80900902?h=60&w=60&fit=cover&mask=circle"></img></a>
<a href="https://github.com/utmcontent"><img src="https://images.weserv.nl/?url=avatars.githubusercontent.com/u/4705133?h=60&w=60&fit=cover&mask=circle"></img></a>
