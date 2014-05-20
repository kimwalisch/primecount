///
/// @file  S1.cpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>

#include <stdint.h>
#include <vector>

namespace primecount {

/// Calculate the contribution of the ordinary leaves in the
/// Lagarias-Miller-Odlyzko algorithm.
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

} // namespace primecount
