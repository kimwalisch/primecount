///
/// @file  pi_lmo4.cpp
/// @brief Implementation of the Lagarias-Miller-Odlyzko prime
///        counting algorithm. This implementation uses the segmented
///        sieve of Eratosthenes and a special tree data structure
///        for faster counting in S2(x).
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
#include <algorithm>
#include <stdint.h>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
  #include "get_omp_threads.hpp"
#endif

using namespace std;

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
  int64_t S2_result = 0;

  std::vector<char> sieve(segment_size);
  std::vector<int32_t> counters(segment_size);
  std::vector<int64_t> next(primes.begin(), primes.end());
  std::vector<int64_t> phi(primes.size(), 0);

  // Segmented sieve of Eratosthenes
  for (int64_t low = 1; low < limit; low += segment_size)
  {
    std::fill(sieve.begin(), sieve.end(), 1);

    // Current segment = interval [low, high[
    int64_t high = std::min(low + segment_size, limit);
    int64_t b = 1;

    // phi(y, b) nodes with b <= c do not contribute to S2, so we
    // simply sieve out the multiples of the first c primes
    for (; b <= c; b++)
    {
      int64_t k = next[b];
      for (int64_t prime = primes[b]; k < high; k += prime)
        sieve[k - low] = 0;
      next[b] = k;
    }

    // Initialize special tree data structure from sieve
    cnt_finit(sieve, counters, segment_size);

    for (; b < pi_y; b++)
    {
      int64_t prime = primes[b];
      int64_t min_m = std::max(x / (prime * high), y / prime);
      int64_t max_m = std::min(x / (prime * low), y);

      // Obviously if (prime >= max_m) then (prime >= lpf[max_m])
      // if so then (prime < lpf[m]) will always evaluate to
      // false and no special leaves are possible
      if (prime >= max_m)
        break;

      for (int64_t m = min_m + 1; m <= max_m; m++)
      {
        if (mu[m] != 0 && prime < lpf[m])
        {
          // We have found a special leaf, compute it's contribution
          // phi(x / (primes[b] * m), b - 1) by counting the number
          // of unsieved elements <= x / (primes[b] * m) after having
          // removed the multiples of the first b - 1 primes
          //
          int64_t n = prime * m;
          int64_t count = cnt_query(counters, (x / n) - low);
          int64_t phi_xn = phi[b] + count;

          S2_result -= mu[m] * phi_xn;
        }
      }

      // Calculate phi(x / ((high - 1) * primes[b]), b) which will
      // be used to calculate special leaves in the next segment
      phi[b] += cnt_query(counters, (high - 1) - low);

      // Remove the multiples of (b)th prime
      int64_t k = next[b];
      for (; k < high; k += prime * 2)
      {
        if (sieve[k - low])
        {
          sieve[k - low] = 0;
          cnt_update(counters, k - low, segment_size);
        }
      }
      next[b] = k;
    }
  }

  return S2_result;
}

} // namespace

namespace primecount {

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo4(int64_t x, int threads)
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

  std::vector<int32_t> lpf = make_least_prime_factor(y);
  std::vector<int32_t> mu = make_moebius(y);
  std::vector<int32_t> primes;
  primes.push_back(0);
  primesieve::generate_primes(y, &primes);

  int64_t pi_y = primes.size() - 1;
  int64_t c = min(PhiTiny::MAX_A, pi_y);
  int64_t s1, s2, p2;

#ifdef _OPENMP
  #pragma omp parallel sections num_threads(get_omp_threads(threads))
  {
    #pragma omp section
    s1 = S1(x, y, c, primes, lpf , mu);
    #pragma omp section
    s2 = S2(x, y, pi_y, c, primes, lpf , mu);
    #pragma omp section
    p2 = P2(x, pi_y, y);
  }
#else
  s1 = S1(x, y, c, primes, lpf , mu);
  s2 = S2(x, y, pi_y, c, primes, lpf , mu);
  p2 = P2(x, pi_y, y);
#endif

  int64_t phi = s1 + s2;
  int64_t sum = phi + pi_y - 1 - p2;

  return sum;
}

} // namespace primecount
