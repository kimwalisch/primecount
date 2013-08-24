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
#include "phi.h"

#include <primesieve/soe/PrimeSieve.h>
#include <primecount.h>
#include <stdint.h>
#include <vector>
#include <limits>

namespace primecount {

/// Initialize a vector with Möbius function values.
/// This implementation is based on code by Rick Sladkey
/// posted here: http://mathoverflow.net/a/99545
///
void init_moebius(std::vector<int32_t>& mu, int64_t max)
{
  mu.resize(max + 1, 1);

  for (int32_t i = 2; i * i <= max; i++)
  {
    if (mu[i] == 1)
    {
      for (int32_t j = i; j <= max; j += i)
        mu[j] *= -i;
      for (int32_t j = i * i; j <= max; j += i * i)
        mu[j] = 0;
    }
  }
  for (int32_t i = 2; i <= max; i++)
  {
    if (mu[i] == i)
      mu[i] = 1;
    else if (mu[i] == -i)
      mu[i] = -1;
    else if (mu[i] < 0)
      mu[i] = 1;
    else if (mu[i] > 0)
      mu[i] = -1;
  }
}

/// Initialize a vector with the least prime
/// factors of the integers <= max.
///
void init_least_factor(std::vector<int32_t>& lpf, int64_t max)
{
  lpf.resize(max + 1, 1);

  // phi(x / 1, c) contributes to the sum,
  // lpf[1] = MAX in order to pass
  // if (least_factor[1] > primes[c])
  if (lpf.size() > 1)
    lpf[1] = std::numeric_limits<int32_t>::max();

  for (int32_t i = 2; i * i <= max; i++)
    if (lpf[i] == 1)
      for (int32_t j = i * 2; j <= max; j += i)
        if (lpf[j] == 1)
          lpf[j] = i;

  for (int32_t i = 2; i <= max; i++)
    if (lpf[i] == 1)
      lpf[i] = i;
}

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3)) operations, O(x^(1/3)) space.
///
int64_t pi_lmo_simple(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  int64_t x13 = iroot<3>(x); 
  int64_t a = pi_meissel(x13);

  // generate the primes <= x^(1/3)
  std::vector<int32_t> primes;
  primes.push_back(0);
  PrimeSieve ps;
  ps.generate_N_Primes(a, &primes);

  std::vector<int32_t> moebius;
  std::vector<int32_t> least_factor;
  init_moebius(moebius, x13);
  init_least_factor(least_factor, x13);

  int64_t c = (a < 7) ? a : 7; 
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
