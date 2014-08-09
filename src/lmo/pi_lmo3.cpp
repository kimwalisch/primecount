///
/// @file  pi_lmo3.cpp
/// @brief Simple implementation of the Lagarias-Miller-Odlyzko prime
///        counting algorithm. This implementation uses the segmented
///        sieve of Eratosthenes to calculate S2(x).
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <pmath.hpp>
#include <generate.hpp>
#include <PhiTiny.hpp>
#include <S1.hpp>

#include <stdint.h>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;

namespace {

/// Calculate the contribution of the special leaves.
/// This implementation uses segmentation which reduces the
/// algorithm's space complexity to O(x^(1/3) * log^2 x).
/// @pre y > 0 && c > 1
///
int64_t S2(int64_t x,
           int64_t y,
           int64_t c,
           vector<int32_t>& primes,
           vector<int32_t>& lpf,
           vector<int32_t>& mu)
{
  int64_t limit = x / y + 1;
  int64_t segment_size = isqrt(limit);
  int64_t pi_y = pi_bsearch(primes, y);
  int64_t S2_result = 0;

  vector<char> sieve(segment_size);
  vector<int64_t> next(primes.begin(), primes.end());
  vector<int64_t> phi(primes.size(), 0);

  // Segmented sieve of Eratosthenes
  for (int64_t low = 1; low < limit; low += segment_size)
  {
    fill(sieve.begin(), sieve.end(), 1);

    // Current segment = interval [low, high[
    int64_t high = min(low + segment_size, limit);

    // phi(y, b) nodes with b <= c do not contribute to S2, so we
    // simply sieve out the multiples of the first c primes
    for (int64_t b = 1; b <= c; b++)
    {
      int64_t k = next[b];
      for (int64_t prime = primes[b]; k < high; k += prime)
        sieve[k - low] = 0;
      next[b] = k;
    }

    for (int64_t b = c + 1; b < pi_y; b++)
    {
      int64_t prime = primes[b];
      int64_t min_m = max(x / (prime * high), y / prime);
      int64_t max_m = min(x / (prime * low), y);
      int64_t i = 0;

      // Obviously if (prime >= max_m) then (prime >= lpf[max_m])
      // if so then (prime < lpf[m]) will always evaluate to
      // false and no special leaves are possible
      if (prime >= max_m)
        break;

      for (int64_t m = max_m; m > min_m; m--)
      {
        if (mu[m] != 0 && prime < lpf[m])
          {
          // We have found a special leaf, compute it's contribution
          // phi(x / (primes[b] * m), b - 1) by counting the number
          // of unsieved elements <= x / (primes[b] * m) after having
          // removed the multiples of the first b - 1 primes
          //
          for (int64_t xn = x / (prime * m); i <= xn - low; i++)
            phi[b] += sieve[i];

          S2_result -= mu[m] * phi[b];
        }
      }

      // Count the remaining unsieved elements in this segment,
      // we need their count in the next segment
      for (; i < high - low; i++)
        phi[b] += sieve[i];

      // Remove the multiples of (b)th prime
      int64_t k = next[b];
      for (; k < high; k += prime * 2)
        sieve[k - low] = 0;
      next[b] = k;
    }
  }

  return S2_result;
}

/// alpha is a tuning factor which should grow like (log(x))^2
/// for the Lagarias-Miller-Odlyzko algorithm
///
double compute_alpha(int64_t x)
{
  double d = (double) x;
  return in_between(1, log(d) * log(d) / 300, iroot<6>(x));
}

} // namespace

namespace primecount {

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * (log x)^2) space.
///
int64_t pi_lmo3(int64_t x)
{
  if (x < 2)
    return 0;

  double alpha = compute_alpha(x);
  int64_t x13 = iroot<3>(x);
  int64_t y = (int64_t) (x13 * alpha);
  int64_t p2 = P2(x, y, 1);

  vector<int32_t> mu = generate_moebius(y);
  vector<int32_t> lpf = generate_least_prime_factors(y);
  vector<int32_t> primes = generate_primes(y);

  int64_t pi_y = primes.size() - 1;
  int64_t c = min(pi_y, PhiTiny::max_a());
  int64_t s1 = S1(x, y, c, primes[c], lpf , mu, 1);
  int64_t s2 = S2(x, y, c, primes, lpf , mu);
  int64_t phi = s1 + s2;
  int64_t sum = phi + pi_y - 1 - p2;

  return sum;
}

} // namespace primecount
