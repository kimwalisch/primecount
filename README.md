primecount
==========
[![Build Status](https://travis-ci.org/kimwalisch/primecount.svg)](https://travis-ci.org/kimwalisch/primecount)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/kimwalisch/primecount?branch=master&svg=true)](https://ci.appveyor.com/project/kimwalisch/primecount)
[![GitHub license](https://img.shields.io/badge/license-BSD%202-blue.svg)](https://github.com/kimwalisch/primecount/blob/master/COPYING)

primecount is a command-line program and C++ library that counts the
primes below an integer x&nbsp;≤&nbsp;10<sup>31</sup> using **highly
optimized** implementations of the
[prime counting function](http://en.wikipedia.org/wiki/Prime-counting_function)
(combinatorial methods). primecount includes implementations of the algorithms
of Legendre, Meissel, Lehmer, Lagarias-Miller-Odlyzko and Deleglise-Rivat all
of which have been parallelized using
[OpenMP](http://en.wikipedia.org/wiki/OpenMP). The Deleglise-Rivat
implementation has 
[recently been distributed](https://github.com/kimwalisch/primecount/blob/master/doc/primecount-MPI.md#primecount-mpi)
using MPI.

primecount contains the **first ever** parallel open source
implementation of the Deleglise-Rivat algorithm and it features a
[novel load balancer](https://github.com/kimwalisch/primecount/blob/master/src/S2LoadBalancer.cpp)
which scales up to hundreds of CPU cores. primecount has already been
used to compute several world records e.g.
[pi(10<sup>27</sup>)](http://www.mersenneforum.org/showthread.php?t=20473) and
[nth_prime(10<sup>24</sup>)](https://oeis.org/A006988), more will hopefully follow!

Binaries
--------
Below are the latest precompiled binaries for Windows 64-bit and Linux x64.
These binaries are statically linked and require a CPU which supports the POPCNT
instruction (2010 or later).

* <a href="https://dl.bintray.com/kimwalisch/primecount/primecount-2.6-win64.zip">primecount-2.6-win64.zip</a>, 386K
* <a href="https://dl.bintray.com/kimwalisch/primecount/primecount-2.6-linux-x64.tar.gz">primecount-2.6-linux-x64.tar.gz</a>, 893K
* Binaries with backup functionality are available [here](https://github.com/kimwalisch/primecount/tree/backup#primecount-backup)

Usage examples
--------------
Open a terminal and run primecount using e.g.:
```sh
# Count the primes below 10^14
$ ./primecount 1e14

# Print progress and status information during computation
$ ./primecount 1e20 --status

# Count primes using Meissel's algorithm
$ ./primecount 2**32 --meissel

# Find the 10^14th prime using 4 threads
$ ./primecount 1e14 --nthprime --threads=4 --time
```

Command-line options
--------------------
```
Usage: primecount x [OPTION]...
Count the primes below x <= 10^31 using fast implementations of the
combinatorial prime counting function.

Options:

  -d,    --deleglise_rivat  Count primes using Deleglise-Rivat algorithm
         --legendre         Count primes using Legendre's formula
         --lehmer           Count primes using Lehmer's formula
  -l,    --lmo              Count primes using Lagarias-Miller-Odlyzko
  -m,    --meissel          Count primes using Meissel's formula
         --Li               Approximate pi(x) using the logarithmic integral
         --Li_inverse       Approximate the nth prime using Li^-1(x)
  -n,    --nthprime         Calculate the nth prime
  -p,    --primesieve       Count primes using the sieve of Eratosthenes
  -s[N], --status[=N]       Show computation progress 1%, 2%, 3%, ...
                            [N] digits after decimal point e.g. N=1, 99.9%
         --test             Run various correctness tests and exit
         --time             Print the time elapsed in seconds
  -t<N>, --threads=<N>      Set the number of threads, 1 <= N <= CPU cores
  -v,    --version          Print version and license information
  -h,    --help             Print this help menu

Advanced Deleglise-Rivat options:

  -a<N>, --alpha=<N>        Tuning factor, 1 <= alpha <= x^(1/6)
         --P2               Only compute the 2nd partial sieve function
         --S1               Only compute the ordinary leaves
         --S2_trivial       Only compute the trivial special leaves
         --S2_easy          Only compute the easy special leaves
         --S2_hard          Only compute the hard special leaves
```

Algorithms
----------
<table>
  <tr>
    <td>Legendre's Formula</td>
    <td><img src="http://kimwalisch.github.io/primecount/formulas/pi_legendre.svg" height="20" align="absmiddle"/></td>
  </tr>
  <tr>
    <td>Meissel's Formula</td>
    <td><img src="http://kimwalisch.github.io/primecount/formulas/pi_meissel.svg" height="20" align="absmiddle"/></td>
  </tr>
  <tr>
    <td>Lehmer's Formula</td>
    <td><img src="http://kimwalisch.github.io/primecount/formulas/pi_lehmer.svg" height="20" align="absmiddle"/></td>
  </tr>
  <tr>
    <td>LMO Formula</td>
    <td><img src="http://kimwalisch.github.io/primecount/formulas/pi_lmo.svg" height="20" align="absmiddle"/></td>
  </tr>
</table>

<p>Up until the early 19th century the most efficient known method for counting
primes was the sieve of Eratosthenes which has a running time of
<img src="http://kimwalisch.github.io/primecount/formulas/Oxloglogx.svg" height="20" align="absmiddle"/>
operations. The first improvement to this bound was Legendre's formula (1830)
which uses the inclusion-exclusion principle to calculate the number of primes
below x without enumerating the individual primes. Legendre's formula has a
running time of
<img src="http://kimwalisch.github.io/primecount/formulas/Ox.svg" height="20" align="absmiddle"/>
operations and uses
<img src="http://kimwalisch.github.io/primecount/formulas/Osqrtx.svg" height="20" align="absmiddle"/>
space. In 1870 E. D. F. Meissel improved Legendre's formula by setting
<img src="http://kimwalisch.github.io/primecount/formulas/apisqrt3x.svg" height="20" align="absmiddle"/>
and by adding the correction term
<img src="http://kimwalisch.github.io/primecount/formulas/P2xa.svg" height="20" align="absmiddle"/>.
Meissel's formula has a running time of
<img src="http://kimwalisch.github.io/primecount/formulas/Omeissel.svg" height="20" align="absmiddle"/>
operations and uses
<img src="http://kimwalisch.github.io/primecount/formulas/Osqrtxlogx.svg" height="20" align="absmiddle"/>
space. In 1959 D. H. Lehmer extended Meissel's formula and slightly improved the running time to
<img src="http://kimwalisch.github.io/primecount/formulas/Olehmer.svg" height="20" align="absmiddle"/>
operations and
<img src="http://kimwalisch.github.io/primecount/formulas/Osqrtxlogx.svg" height="20" align="absmiddle"/>
space. In 1985 J. C. Lagarias, V. S. Miller and A. M. Odlyzko published a new
algorithm based on Meissel's formula which has a lower runtime complexity of
<img src="http://kimwalisch.github.io/primecount/formulas/Oroot23xlogx.svg" height="20" align="absmiddle"/>
operations and which uses only
<img src="http://kimwalisch.github.io/primecount/formulas/Osqrt3xlog2x.svg" height="20" align="absmiddle"/>
space.</p>
<p>For more information on Legendre's, Meissel's and Lehmer's formulas Hans
Riesel's book [4] is probably the best source of information. For the
Lagarias-Miller-Odlyzko algorithm I recommend reading their original paper
[3] as well as Tomás Oliveira's paper [7].

Fast nth prime calculation
--------------------------
The most efficient known method for calculating the nth prime is a combination
of the prime counting function and a prime sieve. The idea is to closely
approximate the nth prime using e.g. the inverse logarithmic integral
<img src="http://kimwalisch.github.io/primecount/formulas/Li-1n.svg" height="20" align="absmiddle"/>
and then count the primes up to this guess using the prime counting function.
Once this is done one starts sieving (e.g. using the segmented sieve of
Eratosthenes) from there on until one finds the actual nth prime. The author
has implemented ```primecount::nth_prime(n)``` this way, it finds the nth
prime in
<img src="http://kimwalisch.github.io/primecount/formulas/Oroot23xlog2x.svg" height="20" align="absmiddle"/>
operations using
<img src="http://kimwalisch.github.io/primecount/formulas/Opisqrtx.svg" height="20" align="absmiddle"/>
space.

Benchmarks
----------
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
  </tr>
  <tr align="right">
    <td>10<sup>11</sup></td>
    <td>4,118,054,813</td>
    <td>0.04s</td>
    <td>0.05s</td>
    <td>0.02s</td>
    <td>0.02s</td>
  </tr>
  </tr>
  <tr align="right">
    <td>10<sup>12</sup></td>
    <td>37,607,912,018</td>
    <td>0.15s</td>
    <td>0.13s</td>
    <td>0.04s</td>
    <td>0.03s</td>
  </tr>
  </tr>
  <tr align="right">
    <td>10<sup>13</sup></td>
    <td>346,065,536,839</td>
    <td>0.70s</td>
    <td>0.46s</td>
    <td>0.11s</td>
    <td>0.07s</td>
  </tr>
  </tr>
  <tr align="right">
    <td>10<sup>14</sup></td>
    <td>3,204,941,750,802</td>
    <td>4.01s</td>
    <td>1.96s</td>
    <td>0.46s</td>
    <td>0.26s</td>
  </tr>
  <tr align="right">
    <td>10<sup>15</sup></td>
    <td>29,844,570,422,669</td>
    <td>27.75s</td>
    <td>12.08s</td>
    <td>1.89s</td>
    <td>1.00s</td>
  </tr>
  <tr align="right">
    <td>10<sup>16</sup></td>
    <td>279,238,341,033,925</td>
    <td>232.30s</td>
    <td>92.09s</td>
    <td>8.20s</td>
    <td>3.77s</td>
  </tr>
  <tr align="right">
    <td>10<sup>17</sup></td>
    <td>2,623,557,157,654,233</td>
    <td>1,836.73s</td>
    <td>731.35s</td>
    <td>36.05s</td>
    <td>14.87s</td>
  </tr>
  <tr align="right">
  <td>10<sup>18</sup></td>
  <td>24,739,954,287,740,860</td>
    <td>14,949.16s</td>
    <td>6,631.73s</td>
    <td>168.46s</td>
    <td>62.03s</td>
  </tr>
  <tr align="right">
  <td>10<sup>19</sup></td>
  <td>234,057,667,276,344,607</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>278.02s</td>
  </tr>
  <tr align="right">
  <td>10<sup>20</sup></td>
  <td>2,220,819,602,560,918,840</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>1,214.14s</td>
  </tr>
  <tr align="right">
  <td>10<sup>21</sup></td>
  <td>21,127,269,486,018,731,928</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>5,422.27s</td>
  </tr>
  <tr align="right">
  <td>10<sup>22</sup></td>
  <td>201,467,286,689,315,906,290</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>NaN</td>
    <td>24,098.16s</td>
  </tr>
</table>

The benchmarks above were run on an Intel Core i7-6700 CPU (4 x 3.4 GHz) from
2015 using a Linux x64 operating system and primecount was compiled using
GCC 5.2.

Build instructions (Unix-like OSes)
-----------------------------------
You need to have installed a C++ compiler and GNU make to build ```primecount```.

```sh
$ ./build.sh
$ sudo make install
```

To build ```primecount``` using
[MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface)
support for distributing computations onto cluster nodes use:
```sh
$ ./build.sh --enable-mpi
```

[primecount-MPI.md](doc/primecount-MPI.md) contains more information.

Build instructions (Microsoft Visual C++)
-----------------------------------------
Open a Visual Studio Command Prompt and run:
```sh
> nmake -f Makefile.msvc
```

C++ API
-------
Below are the main functions declared in
[primesieve.hpp](https://github.com/kimwalisch/primecount/blob/master/include/primecount.hpp).
All functions are multi-threaded by default.

```C++
/// Count the primes <= x
int64_t primecount::pi(int64_t x);

/// 128-bit prime counting function.
/// @param expr  Integer arithmetic expression e.g. "1000", "10^22"
/// @pre   expr  <= 10^31 on 64-bit systems
///        expr    < 2^63 on 32-bit systems
std::string primecount::pi(const std::string& expr);

/// Find the nth prime
int64_t primecount::nth_prime(int64_t n);

int  primecount::get_num_threads();
void primecount::set_num_threads(int threads);
```

C++ library usage
-----------------
Below is an example program that counts the primes below 1000.

```C++
#include <primecount.hpp>
#include <iostream>

int main()
{
    int64_t prime_count = primecount::pi(1000);
    std::cout << "primes below 1000 = " << prime_count << std::endl;
  
    return 0;
}
```

On Unix-like OSes compile using:
```sh
$ c++ -O2 primes.cpp -lprimecount -lprimesieve
```

References
----------
1. A. M. Legendre, Théorie des nombres, Third edition, Paris, 1830. Vol. 2, p. 65.
2. D. H. Lehmer, On the exact number of primes less than a given limit, Illinois J. Math. 3 (1959), pp. 381–388.
3. J. C. Lagarias, V. S. Miller, and A. M. Odlyzko, Computing pi(x): The Meissel-Lehmer method, Mathematics of Computation, 44 (1985), pp. 537–560.
4. M. Deleglise and J. Rivat, "Computing pi(x): The Meissel, Lehmer, Lagarias, Miller, Odlyzko Method", Mathematics of Computation, Volume 65, Number 213, 1996, pp 235–245.
5. Hans Riesel, Prime Numbers and Computer Methods for Factorization, 2nd ed., Birkhäuser, Boston, 1994. pp. 10-38.
6. Raymond Séroul, Programming for Mathematicians, Springer-Verlag, Berlin (2000), pp. 175-181.
7. R. Crandall and C. Pomerance, Prime numbers: a computational perspective, 2nd ed., Springer, New York, 2005. pp. 152-162.
8. Tomás Oliveira e Silva, Computing pi(x): the combinatorial method, Revista do DETUA, vol. 4, no. 6, March 2006, pp. 759-768.
9. Douglas B. Staple, The combinatorial algorithm for computing pi(x), Master of Science Thesis, Dalhousie University Halifax, Nova Scotia, August 2015.
