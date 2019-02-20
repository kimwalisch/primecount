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
  "Usage: primecount x [OPTION]...\n"
  "Count the primes below x <= 10^31 using fast implementations of the\n"
  "combinatorial prime counting function.\n"
  "\n"
  "Backup options:\n"
  "\n"
  "  -b, --backup=<filename>   Set the backup filename. The default backup\n"
  "                            filename is primecount.backup\n"
  "\n"
  "  -r, --resume[=filename]   Resume the last computation from the\n"
  "                            primecount.backup file. If another backup\n"
  "                            filename is provided the computation is resumed\n"
  "                            from that backup file\n"
  "Options:\n"
  "\n"
  "  -d,    --deleglise_rivat  Count primes using Deleglise-Rivat algorithm\n"
  "         --legendre         Count primes using Legendre's formula\n"
  "         --lehmer           Count primes using Lehmer's formula\n"
  "  -m,    --meissel          Count primes using Meissel's formula\n"
  "         --Li               Approximate pi(x) using the logarithmic integral\n"
  "         --Li_inverse       Approximate the nth prime using Li^-1(x)\n"
  "  -n,    --nthprime         Calculate the nth prime\n"
  "  -p,    --primesieve       Count primes using the sieve of Eratosthenes\n"
  "         --phi=<a>          phi(x, a) counts the numbers <= x that are\n"
  "                            not divisible by any of the first a primes\n"
  "         --Ri               Approximate pi(x) using Riemann R\n"
  "         --Ri_inverse       Approximate the nth prime using Ri^-1(x)\n"
  "  -s[N], --status[=N]       Show computation progress 1%, 2%, 3%, ...\n"
  "                            [N] digits after decimal point e.g. N=1, 99.9%\n"
  "         --test             Run various correctness tests and exit\n"
  "         --time             Print the time elapsed in seconds\n"
  "  -t<N>, --threads=<N>      Set the number of threads, 1 <= N <= CPU cores\n"
  "  -v,    --version          Print version and license information\n"
  "  -h,    --help             Print this help menu\n"
  "\n"
  "Advanced Deleglise-Rivat options:\n"
  "\n"
  "  -a<N>, --alpha=<N>        Tuning factor, 1 <= alpha <= x^(1/6)\n"
  "         --P2               Only compute the 2nd partial sieve function\n"
  "         --S1               Only compute the ordinary leaves\n"
  "         --S2_trivial       Only compute the trivial special leaves\n"
  "         --S2_easy          Only compute the easy special leaves\n"
  "         --S2_hard          Only compute the hard special leaves\n"
  "\n"
  "Examples:\n"
  "\n"
  "  primecount 1e13\n"
  "  primecount 1e13 --nthprime --threads=4"
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
