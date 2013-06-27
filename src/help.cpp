///
/// @file  help.cpp
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.h>

#include <iostream>
#include <cstdlib>
#include <string>

using namespace std;

namespace {

const string helpMenu(
  "Usage: primecount x [OPTION]...\n"
  "Count the primes below x < 2^63 using an efficient implementation of the prime\n"
  "counting function e.g. Lehmer's (default) or Meissel's formula.\n"
  "\n"
  "Options:\n"
  "  -l,    --lehmer        Count primes using Lehmer's formula\n"
  "  -g,    --legendre      Count primes using Legendre's formula\n"
  "         --Li            Approximate pi(x) using the logarithmic integral\n"
  "         --Li_inverse    Approximate the nth prime using Li^-1(x)\n"
  "  -m,    --meissel       Count primes using Meissel's formula\n"
  "  -n,    --nth_prime     Calculate the nth prime\n"
  "  -p,    --primesieve    Count primes using the sieve of Eratosthenes\n"
  "         --test          Run various correctness tests and exit\n"
  "  -t<N>, --threads=<N>   Set the number of threads, 1 <= N <= CPU cores\n"
  "  -v,    --version       Print version and license information\n"
  "  -h,    --help          Print this help menu\n"
  "\n"
  "Examples:\n"
  "  primecount 10**13\n"
  "  primecount 10**13 --nth_prime"
);

const string versionInfo(
  "primecount " PRIMECOUNT_VERSION ", <https://github.com/kimwalisch/primecount>\n"
  "Copyright (C) " PRIMECOUNT_YEAR " Kim Walisch\n"
  "BSD 2-Clause License <http://opensource.org/licenses/BSD-2-Clause>"
);

} // end namespace

namespace primecount {

void help()
{
  cout << helpMenu << endl;
  exit(1);
}

void version()
{
  cout << versionInfo << endl;
  exit(1);
}

} // namespace primecount
