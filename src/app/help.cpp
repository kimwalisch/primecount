///
/// @file  help.cpp
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>

#include <iostream>
#include <cstdlib>
#include <string>

using namespace std;

namespace {

const string helpMenu
(
  "Usage: primecount x [options]\n"
  "Count the number of primes less than or equal to x (<= 10^31).\n"
  "\n"
  "Options:\n"
  "\n"
  "  -d,    --deleglise-rivat  Count primes using the Deleglise-Rivat algorithm\n"
  "  -g,    --gourdon          Count primes using Xavier Gourdon's algorithm.\n"
  "                            This is the default algorithm.\n"
  "  -l,    --legendre         Count primes using Legendre's formula\n"
  "         --lehmer           Count primes using Lehmer's formula\n"
  "         --lmo              Count primes using Lagarias-Miller-Odlyzko\n"
  "  -m,    --meissel          Count primes using Meissel's formula\n"
  "         --Li               Approximate pi(x) using the logarithmic integral\n"
  "         --Li-inverse       Approximate the nth prime using Li^-1(x)\n"
  "  -n,    --nth-prime        Calculate the nth prime\n"
  "  -p,    --primesieve       Count primes using the sieve of Eratosthenes\n"
  "         --phi=<a>          phi(x, a) counts the numbers <= x that are not\n"
  "                            divisible by any of the first a primes\n"
  "         --Ri               Approximate pi(x) using Riemann R\n"
  "         --Ri-inverse       Approximate the nth prime using Ri^-1(x)\n"
  "  -s[N], --status[=N]       Show computation progress 1%, 2%, 3%, ...\n"
  "                            Set digits after decimal point: -s1 prints 99.9%\n"
  "         --test             Run various correctness tests and exit\n"
  "         --time             Print the time elapsed in seconds\n"
  "  -t<N>, --threads=<N>      Set the number of threads, 1 <= N <= CPU cores.\n"
  "                            By default primecount uses all available CPU cores.\n"
  "  -v,    --version          Print version and license information\n"
  "  -h,    --help             Print this help menu\n"
  "\n"
  "Advanced options for the Deleglise-Rivat algorithm:\n"
  "\n"
  "  -a<N>, --alpha=<N>        Tuning factor, 1 <= alpha <= x^(1/6)\n"
  "         --P2               Compute the 2nd partial sieve function\n"
  "         --S1               Compute the ordinary leaves\n"
  "         --S2-trivial       Compute the trivial special leaves\n"
  "         --S2-easy          Compute the easy special leaves\n"
  "         --S2-hard          Compute the hard special leaves\n"
  "\n"
  "Advanced options for Xavier Gourdon's algorithm:\n"
  "\n"
  "         --alpha-y=<N>      Tuning factor, 1 <= alpha_y <= x^(1/6)\n"
  "         --alpha-z=<N>      Tuning factor, 1 <= alpha_z <= x^(1/6)\n"
  "         --AC               Compute Gourdon's A + C formulas\n"
  "         --B                Compute Gourdon's B formula\n"
  "         --D                Compute Gourdon's D formula\n"
  "         --Phi0             Compute the Phi0 formula\n"
  "         --Sigma            Compute the 7 Sigma formulas\n"
);

const string versionInfo
(
  "primecount " PRIMECOUNT_VERSION ", <https://github.com/kimwalisch/primecount>\n"
  "Copyright (C) 2013 - 2019 Kim Walisch\n"
  "BSD 2-Clause License <https://opensource.org/licenses/BSD-2-Clause>"
);

} // namespace

namespace primecount {

void help()
{
  cout << helpMenu << endl;
  exit(0);
}

void version()
{
  cout << versionInfo << endl;
  exit(0);
}

} // namespace
