primecount
==========
[![Build Status](https://travis-ci.org/kimwalisch/primecount.svg)](https://travis-ci.org/kimwalisch/primecount)

primecount is a command-line program and C++ library that counts the primes
below an integer x&nbsp;<&nbsp;2<sup>63</sup> using fast implementations
of the prime counting function pi(x). So far primecount offers the option
to count primes using Legendre's, Meissel's and Lehmer's formulas and using
the Lagarias-Miller-Odlyzko algorithm. All implementations have been
parallelized using OpenMP.

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
<img src="http://kimwalisch.github.io/primecount/formulas/Oxlnlnx.svg" height="20" align="absmiddle"/>
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
<img src="http://kimwalisch.github.io/primecount/formulas/Osqrtxlnx.svg" height="20" align="absmiddle"/>
space. In 1959 D. H. Lehmer extended Meissel's formula and slightly improved the running time to
<img src="http://kimwalisch.github.io/primecount/formulas/Olehmer.svg" height="20" align="absmiddle"/>
operations and
<img src="http://kimwalisch.github.io/primecount/formulas/Osqrtxlnx.svg" height="20" align="absmiddle"/>
space. In 1985 J. C. Lagarias, V. S. Miller and A. M. Odlyzko published a new
algorithm based on Meissel's formula which has a lower runtime complexity of
<img src="http://kimwalisch.github.io/primecount/formulas/Oroot23x.svg" height="20" align="absmiddle"/>
operations and which uses only
<img src="http://kimwalisch.github.io/primecount/formulas/Oroot13x.svg" height="20" align="absmiddle"/>
space.</p>
<p>For more information on Legendre's, Meissel's and Lehmer's formulas Hans
Riesel's book [4] is probably the best source of information, for the
Lagarias-Miller-Odlyzko algorithm I recommend reading their original paper
[3] as well as Tomás Oliveira's paper [7].

### Timings

<table>
  <tr align="center">
    <td><b>x</b></td>
    <td><b>Prime Count</b></td>
    <td><b>pi_legendre(x)</b></td>
    <td><b>pi_meissel(x)</b></td>
    <td><b>pi_lehmer(x)</b></td>
    <td><b>pi_lmo(x)</b></td>
  </tr>
  <tr align="right">
    <td>10<sup>10</sup></td>
    <td>455,052,511</td>
    <td>0.05s</td>
    <td>0.04s</td>
    <td>0.03s</td>
    <td>0.02s</td>
  </tr>
  <tr align="right">
    <td>10<sup>11</sup></td>
    <td>4,118,054,813</td>
    <td>0.08s</td>
    <td>0.06s</td>
    <td>0.06s</td>
    <td>0.03s</td>
  </tr>
  <tr align="right">
    <td>10<sup>12</sup></td>
    <td>37,607,912,018</td>
    <td>0.31s</td>
    <td>0.19s</td>
    <td>0.23s</td>
    <td>0.09s</td>
  </tr>
  <tr align="right">
    <td>10<sup>13</sup></td>
    <td>346,065,536,839</td>
    <td>1.49s</td>
    <td>0.95s</td>
    <td>1.03s</td>
    <td>0.33s</td>
  </tr>
  <tr align="right">
    <td>10<sup>14</sup></td>
    <td>3,204,941,750,802</td>
    <td>9.30s</td>
    <td>6.81s</td>
    <td>5.05s</td>
    <td>1.26s</td>
  </tr>
  <tr align="right">
    <td>10<sup>15</sup></td>
    <td>29,844,570,422,669</td>
    <td>60.38s</td>
    <td>49.12s</td>
    <td>28.26s</td>
    <td>5.94s</td>
  </tr>
  <tr align="right">
    <td>10<sup>16</sup></td>
    <td>279,238,341,033,925</td>
    <td>423.12s</td>
    <td>337.56s</td>
    <td>173.78s</td>
    <td>29.20s</td>
  </tr>
  <tr align="right">
    <td>10<sup>17</sup></td>
    <td>2,623,557,157,654,233</td>
    <td>3749.72s</td>
    <td>2879.81s</td>
    <td>1335.85s</td>
    <td>156.48s</td>
  </tr>
</table>

The benchmarks above were run on an Intel Core i7-4770 CPU (4 x 3.4GHz) from
2013 using a 64-bit Linux operating system. primecount was compiled using GCC
4.8 and used all CPU cores for each benchmark.

### Fast nth prime calculation

The most efficient known method for calculating the nth prime is a combination
of the prime counting function and a prime sieve. The idea is to closely
approximate the nth prime using e.g. the inverse logarithmic integral
<img src="http://kimwalisch.github.io/primecount/formulas/Li-1n.svg" height="20" align="absmiddle"/>
and then count the primes up to this guess using the prime counting function.
Once this is done one starts sieving (e.g. using the segmented sieve of
Eratosthenes) from there on until one finds the actual nth prime. The author
has implemented ```primecount::nth_prime(n)``` this way. In practice most time
is spent by the prime counting function so the calculation of the nth prime is
about as fast as counting the primes below the nth prime.

### Precompiled binaries

Below are the latest precompiled binaries for Windows 64-bit and Linux x86-64.
These binaries are command-line programs for use with a terminal or Command
Prompt.

* <a href="http://dl.bintray.com/kimwalisch/primecount/primecount-0.18-win64.zip">primecount-0.18-win64.zip</a>,  	214K
* <a href="http://dl.bintray.com/kimwalisch/primecount/primecount-0.18-linux-x64.tar.gz">primecount-0.18-linux-x64.tar.gz</a>, 662K

SHA1 checksums of the files:
```sh
7f316e836bcf99145c256fdd875a79d0ffc5df8c  primecount-0.18-win64.zip
873ae49dc63264190bddd272a7e18260550fd204  primecount-0.18-linux-x64.tar.gz
```

### Usage examples
Open a terminal and run the primecount command-line application using e.g.:
```sh
# Count the primes below 10^14
$ ./primecount 10**14

# Count the primes below 10^14 using Meissel's algorithm
$ ./primecount 10**14 --meissel

# Find the 10^14th prime
$ ./primecount 10**14 --nthprime

# Print an option summary
$ ./primecount --help
```

### Build instructions (Unix-like OSes)
In order to build primecount you need to have installed a C++ compiler and
GNU make (or any other make). primecount depends on the author's primesieve
library, download it from http://primesieve.org/downloads
and install it using:
```sh
$ ./configure
$ make
$ sudo make install
```
Then download the latest
<a href="http://dl.bintray.com/kimwalisch/primecount/primecount-0.18.tar.gz">primecount-0.18.tar.gz</a>
release tarball and build it using:
```sh
$ ./configure
$ make
$ sudo make install
```
If you have installed primesieve but primecount's configure script still fails
due to missing libprimesieve then you need to add /usr/local/lib to your library
path:
```sh
export LIBRARY_PATH=/usr/local/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
```

### Build instructions (Microsoft Visual C++)
In order to facilitate building primecount on Windows the build process has been
completely automated i.e. the ```nmake``` command automatically downloads and
builds the latest primesieve library (dependency) before building primecount. To
build primecount simply open a Visual Studio Command Prompt and run:
```sh
> nmake -f Makefile.msvc
```

### C++ API
Below is a list of the functions declared in the ````primecount.hpp```` header
file. A short description of each function including its run-time and space
complexity can be read <a href="include/primecount.hpp">here</a>.

```C++
/// @file  primecount.hpp

int64_t primecount::pi            (int64_t x, int threads = MAX_THREADS);
int64_t primecount::pi_legendre   (int64_t x, int threads = MAX_THREADS);
int64_t primecount::pi_lehmer     (int64_t x, int threads = MAX_THREADS);
int64_t primecount::pi_lmo        (int64_t x, int threads = MAX_THREADS);
int64_t primecount::pi_meissel    (int64_t x, int threads = MAX_THREADS);
int64_t primecount::pi_primesieve (int64_t x, int threads = MAX_THREADS);
int64_t primecount::nth_prime     (int64_t n, int threads = MAX_THREADS);
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
4. Hans Riesel, Prime Numbers and Computer Methods for Factorization, 2nd ed., Birkhäuser, Boston, 1994. pp. 10-38.
5. Raymond Séroul, Programming for Mathematicians, Springer-Verlag, Berlin (2000), pp. 175-181.
6. R. Crandall and C. Pomerance, Prime numbers: a computational perspective, 2nd ed., Springer, New York, 2005. pp. 152-162.
7. Tomás Oliveira e Silva, Computing pi(x): the combinatorial method, Revista do DETUA, vol. 4, no. 6, March 2006, pp. 759-768.
