primecount
==========
[![Build Status](https://travis-ci.org/kimwalisch/primecount.svg)](https://travis-ci.org/kimwalisch/primecount)

primecount is a command-line program and C++ library that counts the primes
below an integer x&nbsp;≤&nbsp;10<sup>27</sup> using highly optimized
implementations of the prime counting function pi(x). So far primecount
offers the option to count primes using the algorithms of Legendre, Meissel,
Lehmer, Lagarias-Miller-Odlyzko and Deleglise-Rivat. All implementations
have been parallelized using OpenMP.

### Algorithms and complexity

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

### Timings

<table>
  <tr align="center">
    <td><b>x</b></td>
    <td><b>Prime Count</b></td>
    <td><b>Legendre</b></td>
    <td><b>Lehmer</b></td>
    <td><b>Lagarias<br/>Miller<br/>Odlyzko</b></td>
    <td><b>Deleglise<br/>Rivat</b></td>
  </tr>
  <tr align="right">
    <td>10<sup>10</sup></td>
    <td>455,052,511</td>
    <td>0.03s</td>
    <td>0.02s</td>
    <td>0.01s</td>
    <td>0.00s</td>
  </tr>
  </tr>
  <tr align="right">
    <td>10<sup>11</sup></td>
    <td>4,118,054,813</td>
    <td>0.07s</td>
    <td>0.05s</td>
    <td>0.02s</td>
    <td>0.01s</td>
  </tr>
  </tr>
  <tr align="right">
    <td>10<sup>12</sup></td>
    <td>37,607,912,018</td>
    <td>0.28s</td>
    <td>0.21s</td>
    <td>0.07s</td>
    <td>0.04s</td>
  </tr>
  </tr>
  <tr align="right">
    <td>10<sup>13</sup></td>
    <td>346,065,536,839</td>
    <td>1.35s</td>
    <td>0.99s</td>
    <td>0.30s</td>
    <td>0.13s</td>
  </tr>
  </tr>
  <tr align="right">
    <td>10<sup>14</sup></td>
    <td>3,204,941,750,802</td>
    <td>8.52s</td>
    <td>4.89s</td>
    <td>1.20s</td>
    <td>0.47s</td>
  </tr>
  <tr align="right">
    <td>10<sup>15</sup></td>
    <td>29,844,570,422,669</td>
    <td>56.22s</td>
    <td>27.167s</td>
    <td>5.17s</td>
    <td>1.70s</td>
  </tr>
  <tr align="right">
    <td>10<sup>16</sup></td>
    <td>279,238,341,033,925</td>
    <td>387.03s</td>
    <td>170.47s</td>
    <td>23.61s</td>
    <td>6.44s</td>
  </tr>
  <tr align="right">
    <td>10<sup>17</sup></td>
    <td>2,623,557,157,654,233</td>
    <td>3,749.72s</td>
    <td>1,335.85s</td>
    <td>114.02s</td>
    <td>25.54s</td>
  </tr>
  <tr align="right">
  <td>10<sup>18</sup></td>
  <td>24,739,954,287,740,860</td>
    <td>31,897.66s</td>
    <td>9,885.71s</td>
    <td>743.98s</td>
    <td>102.10s</td>
  </tr>
</table>

The benchmarks above were run on an Intel Core i7-4770 CPU (4 x 3.4 GHz) from
2013 using a 64-bit Linux operating system and primecount was compiled using
GCC 4.8.

<table>
  <tr align="center">
    <td><b>x</b></td>
    <td><b>Prime Count</b></td>
    <td><b>Deleglise<br/>Rivat</b></td>
    <td><b>Threads</b></td>
    <td><b>Memory Usage</b></td>
  </tr>
  <tr align="right">
    <td>10<sup>19</sup></td>
    <td>234,057,667,276,344,607</td>
    <td>1.74 min</td>
    <td>88</td>
    <td>0.22 GB</td>
  </tr>
  <tr align="right">
    <td>10<sup>20</sup></td>
    <td>2,220,819,602,560,918,840</td>
    <td>7.15 min</td>
    <td>88</td>
    <td>0.54 GB</td>
  </tr>
  <tr align="right">
    <td>10<sup>21</sup></td>
    <td>21,127,269,486,018,731,928</td>
    <td>39.81 min</td>
    <td>44</td>
    <td>1.11 GB</td>
  </tr>
  <tr align="right">
    <td>10<sup>22</sup></td>
    <td>201,467,286,689,315,906,290</td>
    <td>5.58 hours</td>
    <td>44</td>
    <td>2.63 GB</td>
  </tr>
</table>

The benchmarks above were run in November 2014 on a system with 22 IBM POWER8
CPU cores running at 3.0 GHz with 8 threads per CPU core. The operating system
used was Linux and primecount was compiled using GCC 4.8. Above
3&nbsp;*&nbsp;10<sup>22</sup> primecount's memory usage triples because it must
then use 64-bit arrays (instead of a 32-bit arrays) for some of its data
structures.

### Fast nth prime calculation

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

### Precompiled binaries

Below are the latest precompiled binaries for Windows 64-bit and Linux x86-64.
These binaries are statically linked and require a CPU (2010 or later) which
supports the POPCNT instruction.

* <a href="http://dl.bintray.com/kimwalisch/primecount/primecount-1.2-win64.zip">primecount-1.2-win64.zip</a>, 355K
* <a href="http://dl.bintray.com/kimwalisch/primecount/primecount-1.2-linux-x64.tar.gz">primecount-1.2-linux-x64.tar.gz</a>, 114K

SHA1 checksums of the files:
```sh
5e68982b272f628285d19be75c9adf7f72c4ff27  primecount-1.2-win64.zip
2fe46a048a0da2a33bed933dbd3d75feb76a9285  primecount-1.2-linux-x64.tar.gz
```

### Usage examples
Open a terminal and run the primecount command-line application using e.g.:
```sh
# Count the primes below 10^14
$ ./primecount 1e14

# Print progress and status information during computation
$ ./primecount 1e20 --status

# Count the primes using Meissel's algorithm
$ ./primecount 1e14 --meissel

# Find the 10^14th prime
$ ./primecount 1e14 --nthprime --time

# Print an option summary
$ ./primecount --help
```

### Build instructions (Unix-like OSes)
To build primecount you need to have installed a C++ compiler and GNU make.
primecount depends on the author's primesieve library, download it from
http://primesieve.org/downloads and install it using:
```sh
$ ./configure
$ make
$ sudo make install
```
If you are not using Linux then you may need to export these variables:
```sh
export LIBRARY_PATH=/usr/local/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
export CPLUS_INCLUDE_PATH=/usr/local/include:$CPLUS_INCLUDE_PATH
```
Finally download the latest
<a href="http://dl.bintray.com/kimwalisch/primecount/primecount-1.2.tar.gz">primecount-1.2.tar.gz</a>
release tarball and build it using:
```sh
$ ./configure
$ make
$ sudo make install
```
If your CPU supports the
[POPCNT instruction](http://en.wikipedia.org/wiki/SSE4#POPCNT_and_LZCNT)
then it is enabled in the build process. Using POPCNT speeds up
primecount by about 10 percent. If you need maximum portability you can
disable POPCNT using:
```sh
$ ./configure --disable-popcnt
```

### Build instructions (Microsoft Visual C++)
In order to facilitate building primecount on Windows the build process has been
completely automated i.e. the ```nmake``` command automatically downloads and
builds the latest primesieve library (dependency) before building primecount. To
build primecount simply open a Visual Studio Command Prompt and run:
```sh
> nmake -f Makefile.msvc
```

### 128-bit support
If you have a 64-bit CPU and your compiler supports the ```__int128_t```
type (GCC, Clang, Cygwin64, Mingw-w64) then 128-bit support is automatically
enabled in the build process and you can count primes up to 10<sup>27</sup>.

### C++ API
Below is a list of the functions declared in the ````primecount.hpp```` header
file. A short description of each function including its run-time and space
complexity can be read <a href="include/primecount.hpp">here</a>.

```C++
/// @file  primecount.hpp

int64_t primecount::pi                 (int64_t x);
int64_t primecount::pi_deleglise_rivat (int64_t x);
int64_t primecount::pi_legendre        (int64_t x);
int64_t primecount::pi_lehmer          (int64_t x);
int64_t primecount::pi_lmo             (int64_t x);
int64_t primecount::pi_meissel         (int64_t x);
int64_t primecount::pi_primesieve      (int64_t x);
int64_t primecount::nth_prime          (int64_t n);

/// 128-bit prime counting function.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
/// @param expr  Integer arithmetic expression e.g. "1000", "10^22"
/// @pre   expr  <= primecount::max()
///
std::string primecount::pi(const std::string& expr);

/// @return  Largest integer supported by pi(const std::string&)
std::string primecount::max();

int  primecount::get_num_threads();
void primecount::set_num_threads(int threads);
```

### Using libprimecount
Below is a C++ example program that counts the primes below n and calculates
the nth prime using libprimecount.

```C++
#include <primecount.hpp>
#include <iostream>
#include <cstdlib>

int main(int, char** argv)
{
  int64_t n = 1000;
  if (argv[1])
    n = std::atol(argv[1]);

  int64_t prime_count = primecount::pi(n);
  std::cout << "primes below " << n << " = " << prime_count << std::endl;

  int64_t nth_prime = primecount::nth_prime(n);
  std::cout << n << "th prime = " << nth_prime << std::endl;

  return 0;
}
```

On Unix-like operating systems compile using:
```sh
$ c++ -O2 primes.cpp -lprimecount
```

### References
1. A. M. Legendre, Théorie des nombres, Third edition, Paris, 1830. Vol. 2, p. 65.
2. D. H. Lehmer, On the exact number of primes less than a given limit, Illinois J. Math. 3 (1959), pp. 381–388.
3. J. C. Lagarias, V. S. Miller, and A. M. Odlyzko, Computing pi(x): The Meissel-Lehmer method, Mathematics of Computation, 44 (1985), pp. 537–560.
4. M. Deleglise and J. Rivat, "Computing pi(x): The Meissel, Lehmer, Lagarias, Miller, Odlyzko Method", Mathematics of Computation, Volume 65, Number 213, 1996, pp 235–245.
5. Hans Riesel, Prime Numbers and Computer Methods for Factorization, 2nd ed., Birkhäuser, Boston, 1994. pp. 10-38.
6. Raymond Séroul, Programming for Mathematicians, Springer-Verlag, Berlin (2000), pp. 175-181.
7. R. Crandall and C. Pomerance, Prime numbers: a computational perspective, 2nd ed., Springer, New York, 2005. pp. 152-162.
8. Tomás Oliveira e Silva, Computing pi(x): the combinatorial method, Revista do DETUA, vol. 4, no. 6, March 2006, pp. 759-768.
