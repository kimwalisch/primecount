///
/// @file  pi_lmo5.cpp
/// @brief Work in progress...
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "PhiTiny.hpp"
#include "Pk.hpp"
#include "pmath.hpp"
#include "tos_counters.hpp"

#include <primecount.hpp>
#include <primesieve.hpp>
#include <stdint.h>
#include <vector>
#include <iostream>

using std::log;

namespace {

/// Calculate the contribution of the ordinary leaves.
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
      S1_result += mu[n] * primecount::phi(x / n, c);

  return S1_result;
}

/// Calculate the contribution of the special leaves.
/// @pre c >= 2
///
int64_t S2(int64_t x,
           int64_t y,
           int64_t pi_y,
           int64_t c,
           std::vector<int32_t>& primes,
           std::vector<int32_t>& lpf,
           std::vector<int32_t>& mu)
{
  int64_t limit = x / y + 1;
  int64_t segment_size = next_power_of_2(isqrt(limit));
  int64_t sqrt_y = isqrt(y);
  int64_t S2_result = 0;

  // vector used for sieving
  std::vector<char> sieve(segment_size);
  std::vector<int32_t> counters(segment_size);
  std::vector<int32_t> pi = make_pi(y);
  std::vector<int32_t> l_max(primes.size(), primes.size() - 1);
  std::vector<int64_t> next(primes.begin(), primes.end());
  std::vector<int64_t> phi(primes.size(), 0);

  // Segmented sieve of Eratosthenes
  for (int64_t low = 1; low < limit; low += segment_size)
  {
    std::fill(sieve.begin(), sieve.end(), 1);

    // Current segment = interval [low, high[
    int64_t high = std::min(low + segment_size, limit);

    // phi(y, b) nodes with b <= c do not contribute to S2, so we
    // simply sieve out the multiples of the first c primes
    for (int64_t b = 1; b <= c; b++)
    {
      int64_t k = next[b];
      for (int64_t prime = primes[b]; k < high; k += prime)
        sieve[k - low] = 0;
      next[b] = k;
    }

    // Initialize special tree data structure from sieve
    cnt_finit(sieve, counters, segment_size);

    int64_t b = c;
    for (; b + 1 < pi[sqrt_y]; b++)
    {
      int64_t prime = primes[b + 1];
      int64_t min_m = std::max(x / (high * prime), y / prime);
      int64_t max_m = std::min(x / (low * prime), y);

      if (prime >= max_m)
        break;

      for (int64_t m = max_m; m > min_m; m--)
      {
        if (mu[m] != 0 && prime < lpf[m])
        {
          int64_t n = prime * m;
          int64_t count = cnt_query(counters, (x / n) - low);
          int64_t phi_xn = phi[b + 1] + count;

          S2_result -= mu[m] * phi_xn;
        }
      }

      phi[b + 1] += cnt_query(counters, (high - 1) - low);

      int64_t k = next[b + 1];
      for (; k < high; k += prime * 2)
      {
        if (sieve[k - low])
        {
          sieve[k - low] = 0;
          cnt_update(counters, k - low, segment_size);
        }
      }
      next[b + 1] = k;
    }

    int64_t special_leaf_threshold = std::max(x / high, y);

    for (; b + 1 < pi_y; b++)
    {
      int64_t prime = primes[b + 1];
      int64_t l = l_max[b + 1];
      special_leaf_threshold = std::max(prime * prime, special_leaf_threshold);

      if (prime >= primes[l])
        break;

      for (; prime * primes[l] > special_leaf_threshold; l--)
      {
        int64_t n = prime * primes[l];
        int64_t phi_xn = phi[b + 1] + cnt_query(counters, (x / n) - low);
        S2_result += phi_xn;
      }

      l_max[b + 1] = l;
      phi[b + 1] += cnt_query(counters, (high - 1) - low);

      int64_t k = next[b + 1];
      for (; k < high; k += prime * 2)
      {
        if (sieve[k - low])
        {
          sieve[k - low] = 0;
          cnt_update(counters, k - low, segment_size);
        }
      }
      next[b + 1] = k;
    }
  }

  return S2_result;
}

} // namespace

namespace primecount {

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3)) operations, O(x^0.5) space.
/// @note O(x^0.5) space is due to parallel P2(x, a).
///
int64_t pi_lmo5(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  // Optimization factor, see:
  // J. C. Lagarias, V. S. Miller, and A. M. Odlyzko, Computing pi(x): The Meissel-
  // Lehmer method, Mathematics of Computation, 44 (1985), p. 556.
  double beta = 0.6;
  double alpha = std::max(1.0, log(log((double) x)) * beta);

  int64_t x13 = iroot<3>(x);
  int64_t y = (int64_t)(x13 * alpha);
  int64_t a = pi_lehmer(y);
  int64_t c = (a < 6) ? a : 6;

  std::vector<int32_t> lpf = make_least_prime_factor(y);
  std::vector<int32_t> mu = make_moebius(y);
  std::vector<int32_t> primes;
  primes.push_back(0);
  primesieve::generate_n_primes(a, &primes);

  int64_t phi = S1(x, y, c, primes, lpf , mu) + S2(x, y, a, c, primes, lpf , mu);
  int64_t sum = phi + a - 1 - P2(x, a, threads);

  return sum;
}

} // namespace primecount
