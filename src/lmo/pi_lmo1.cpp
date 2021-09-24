///
/// @file  pi_lmo1.cpp
/// @brief Simple demonstration implementation of the
///        Lagarias-Miller-Odlyzko prime counting algorithm.
///        Usually in the Lagarias-Miller-Odlyzko algorithm phi(x, a)
///        is calculated using a prime sieve but this simple
///        implementation calculates phi(x, a) using the recursive
///        formula with caching.
///
///        Lagarias-Miller-Odlyzko formula:
///        pi(x) = pi(y) + S1(x, a) + S2(x, a) - 1 - P2(x, a)
///        with y = x^(1/3), a = pi(y)
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <generate.hpp>
#include <PhiTiny.hpp>
#include <imath.hpp>

#include <stdint.h>

namespace primecount {

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3))
/// Memory usage: O(x^(1/2))
///
int64_t pi_lmo1(int64_t x)
{
  if (x < 2)
    return 0;

  int threads = 1;
  int64_t y = iroot<3>(x);
  int64_t pi_y = pi_noprint(y, threads);
  int64_t c = PhiTiny::get_c(y);
  int64_t S1 = 0;
  int64_t S2 = 0;

  auto primes = generate_primes<int32_t>(y);
  auto lpf = generate_lpf(y);
  auto mu = generate_moebius(y);
  bool is_print = false;

  // ordinary leaves
  for (int64_t n = 1; n <= y; n++)
    if (lpf[n] > primes[c])
      S1 += mu[n] * phi(x / n, c, threads, is_print);

  // special leaves
  for (int64_t b = c + 1; b < pi_y; b++)
    for (int64_t m = (y / primes[b]) + 1; m <= y; m++)
      if (lpf[m] > primes[b])
        S2 -= mu[m] * phi(x / (primes[b] * m), b - 1, threads, is_print);

  int64_t phi = S1 + S2;
  int64_t sum = phi + pi_y - 1 - P2(x, y, threads);

  return sum;
}

} // namespace
