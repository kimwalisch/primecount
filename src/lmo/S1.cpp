///
/// @file  S1.cpp
/// @brief Functions to calculate the contribution of the ordinary
///        leaves in the Lagarias-Miller-Odlyzko algorithm.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <FactorTable.hpp>
#include <pmath.hpp>

#include <stdint.h>
#include <vector>

namespace primecount {

/// Run time: O(y) operations, O(y) space.
///
int64_t S1(int64_t x,
           int64_t y,
           int64_t c,
           std::vector<int32_t>& primes,
           std::vector<int32_t>& lpf,
           std::vector<int32_t>& mu)
{
  int64_t S1_result = 0;

  for (int64_t n = 1; n <= y; n++)
    if (lpf[n] > primes[c])
      S1_result += mu[n] * phi(x / n, c);

  return S1_result;
}

/// This version uses 17 times less memory than the version above.
/// Run time: O(y) operations, O(y) space.
///
int64_t S1(int64_t x,
           int64_t y,
           int64_t c,
           std::vector<int32_t>& primes,
           FactorTable& factors)
{
  // FactorTable contains only numbers coprime to 2, 3, 5 and 7
  if (primes[c] <= 7)
  {
    std::vector<int32_t> mu = make_moebius(y);
    std::vector<int32_t> lpf = make_least_prime_factor(y);
    return S1(x, y, c, primes, lpf, mu);
  }

  int64_t S1_result = 0;
  int64_t limit = factors.get_index(y);

  for (int64_t i = factors.get_index(1); i <= limit; i++)
    if (factors.lpf(i) > primes[c])
      S1_result += factors.mu(i) * phi(x / factors.get_number(i), c);

  return S1_result;
}

} // namespace primecount
