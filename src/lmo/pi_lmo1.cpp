///
/// @file  pi_lmo1.cpp
/// @brief Simple demonstration implementation of the
///        Lagarias-Miller-Odlyzko prime counting algorithm.
///        Usually in the Lagarias-Miller-Odlyzko algorithm phi(x, a)
///        is calculated using a prime sieve but this simple
///        implementation calculates phi(x, a) using the recursive
///        formula with caching.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <generate.hpp>
#include <pmath.hpp>
#include <PhiCache.hpp>

#include <stdint.h>
#include <algorithm>
#include <vector>

using namespace std;

namespace primecount {

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3)) operations, O(x^(1/3)) space.
///
int64_t pi_lmo1(int64_t x)
{
  if (x < 2)
    return 0;

  int64_t y = iroot<3>(x); 
  int64_t pi_y = pi_meissel(y, 1);
  int64_t c = min(pi_y, PhiTiny::max_a());
  int64_t S1 = 0;
  int64_t S2 = 0;

  vector<int32_t> lpf = generate_least_prime_factors(y);
  vector<int32_t> mu = generate_moebius(y);
  vector<int32_t> primes = generate_primes(y);

  // Calculate the contribution of the ordinary leaves
  for (int64_t n = 1; n <= y; n++)
    if (lpf[n] > primes[c])
      S1 += mu[n] * phi(x / n, c);

  PhiCache cache(primes);

  // Calculate the contribution of the special leaves
  for (int64_t b = c + 1; b < pi_y; b++)
    for (int64_t m = (y / primes[b]) + 1; m <= y; m++)
      if (lpf[m] > primes[b])
        S2 -= mu[m] * phi(x / (primes[b] * m), b - 1, &cache);

  int64_t phi = S1 + S2;
  int64_t sum = phi + pi_y - 1 - P2(x, y, 1);

  return sum;
}

} // namespace primecount
