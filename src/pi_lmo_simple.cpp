///
/// @file  pi_lmo_simple.cpp
/// @brief Simple demonstration implementation of the
///        Lagarias-Miller-Odlyzko prime counting algorithm.
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "imath.h"
#include "Pk.h"

#include <primecount.h>
#include <primesieve.hpp>
#include <stdint.h>
#include <vector>

namespace primecount {

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3)) operations, O(x^0.5) space.
/// @note O(x^0.5) space is due to parallel P2(x, a).
///
int64_t pi_lmo_simple(int64_t x, int)
{
  if (x < 2)
    return 0;

  int64_t x13 = iroot<3>(x); 
  int64_t a = pi_meissel(x13);

  // generate the primes <= x^(1/3)
  std::vector<int32_t> primes;
  primes.push_back(0);
  primesieve::generate_n_primes(a, &primes);

  std::vector<int32_t> moebius;
  std::vector<int32_t> least_factor;
  init_moebius(moebius, x13);
  init_least_factor(least_factor, x13);

  int64_t c = (a < 6) ? a : 6;
  int64_t S1 = 0;
  int64_t S2 = 0;

  for (int64_t n = 1; n <= x13; n++)
    if (least_factor[n] > primes[c])
      S1 += moebius[n] * phi(x / n, c);

  for (int64_t b = c; b + 1 < a; b++)
    for (int64_t m = (x13 / primes[b + 1]) + 1; m <= x13; m++)
      if (least_factor[m] > primes[b + 1])
        S2 -= moebius[m] * phi(x / (m * primes[b + 1]), b);

  int64_t phi = S1 + S2;
  int64_t sum = phi + a - 1 - P2(x, a);

  return sum;
}

} // namespace primecount
