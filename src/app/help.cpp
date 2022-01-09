///
/// @file  help.cpp
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>

#include <iostream>
#include <cstdlib>
#include <string>

namespace primecount {

void help(int exitCode)
{
  const std::string helpMenu =
    "Usage: primecount x [options]\n"
    "Count the number of primes less than or equal to x (<= 10^31).\n"
    "\n"
    "Options:\n"
    "\n"
    "  -d, --deleglise-rivat  Count primes using the Deleglise-Rivat algorithm\n"
    "  -g, --gourdon          Count primes using Xavier Gourdon's algorithm.\n"
    "                         This is the default algorithm.\n"
    "  -l, --legendre         Count primes using Legendre's formula\n"
    "      --lehmer           Count primes using Lehmer's formula\n"
    "      --lmo              Count primes using Lagarias-Miller-Odlyzko\n"
    "  -m, --meissel          Count primes using Meissel's formula\n"
    "      --Li               Approximate pi(x) using the logarithmic integral\n"
    "      --Li-inverse       Approximate the nth prime using Li^-1(x)\n"
    "  -n, --nth-prime        Calculate the nth prime\n"
    "  -p, --primesieve       Count primes using the sieve of Eratosthenes\n"
    "      --phi <X> <A>      phi(x, a) counts the numbers <= x that are not\n"
    "                         divisible by any of the first a primes\n"
    "      --Ri               Approximate pi(x) using Riemann R\n"
    "      --Ri-inverse       Approximate the nth prime using Ri^-1(x)\n"
    "  -s, --status[=NUM]     Show computation progress 1%, 2%, 3%, ...\n"
    "                         Set digits after decimal point: -s1 prints 99.9%\n"
    "      --test             Run various correctness tests and exit\n"
    "      --time             Print the time elapsed in seconds\n"
    "  -t, --threads=NUM      Set the number of threads, 1 <= NUM <= CPU cores.\n"
    "                         By default primecount uses all available CPU cores.\n"
    "  -v, --version          Print version and license information\n"
    "  -h, --help             Print this help menu\n"
    "\n"
    "Advanced options for the Deleglise-Rivat algorithm:\n"
    "\n"
    "  -a, --alpha=NUM        Set tuning factor: y = x^(1/3) * alpha\n"
    "      --P2               Compute the 2nd partial sieve function\n"
    "      --S1               Compute the ordinary leaves\n"
    "      --S2-trivial       Compute the trivial special leaves\n"
    "      --S2-easy          Compute the easy special leaves\n"
    "      --S2-hard          Compute the hard special leaves\n"
    "\n"
    "Advanced options for Xavier Gourdon's algorithm:\n"
    "\n"
    "      --alpha-y=NUM      Set tuning factor: y = x^(1/3) * alpha_y\n"
    "      --alpha-z=NUM      Set tuning factor: z = y * alpha_z\n"
    "      --AC               Compute the A + C formulas\n"
    "      --B                Compute the B formula\n"
    "      --D                Compute the D formula\n"
    "      --Phi0             Compute the Phi0 formula\n"
    "      --Sigma            Compute the 7 Sigma formulas\n";

  std::cout << helpMenu << std::endl;
  std::exit(exitCode);
}

void version()
{
  const std::string versionInfo =
    "primecount " PRIMECOUNT_VERSION ", <https://github.com/kimwalisch/primecount>\n"
    "Copyright (C) 2013 - 2022 Kim Walisch\n"
    "BSD 2-Clause License <https://opensource.org/licenses/BSD-2-Clause>";

  std::cout << versionInfo << std::endl;
  std::exit(0);
}

} // namespace
