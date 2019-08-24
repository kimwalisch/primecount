///
/// @file  D1.cpp
/// @brief Simple demonstration implementation of the D(x, y) formula
///        in Xavier Gourdon's prime counting algorithm. This
///        implementation runs single threaded and does not use the
///        highly optimized Sieve.cpp.
///
///        The D formula corresponds to the computation of the hard
///        special leaves (those that require use of a sieve) in the
///        Lagarias-Miller-Odlyzko and Deleglise-Rivat prime counting
///        algorithms.
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <gourdon.hpp>
#include <primecount-internal.hpp>
#include <generate.hpp>
#include <min.hpp>
#include <imath.hpp>
#include <print.hpp>

#include <stdint.h>
#include <vector>

using namespace std;

namespace primecount {

int64_t D(int64_t x,
          int64_t y,
          int64_t z,
          int64_t k)
{
  print("");
  print("=== D(x, y) ===");
  print_gourdon(x, y, z, k, 1);

  double time = get_time();
  int64_t sum = 0;
  int64_t limit = x / z + 1;
  int64_t segment_size = isqrt(limit);
  int64_t x_star = get_x_star_gourdon(x, y);

  auto pi = generate_pi(x_star);
  auto primes = generate_primes<int32_t>(x_star);
  auto mu = generate_moebius(z);
  auto lpf = generate_lpf(z);
  auto mpf = generate_mpf(z);

  vector<char> sieve(segment_size);
  vector<int64_t> next(primes.begin(), primes.end());
  vector<int64_t> phi(primes.size(), 0);
  int64_t pi_x_star = pi[x_star];

  // segmented sieve of Eratosthenes
  for (int64_t low = 1; low < limit; low += segment_size)
  {
    // Current segment [low, high[
    int64_t high = min(low + segment_size, limit);
    int64_t b = 1;

    // Reset the sieve array
    fill(sieve.begin(), sieve.end(), 1);

    // Pre-sieve multiples of first k primes
    for (; b <= k; b++)
    {
      int64_t j = next[b];
      for (int64_t prime = primes[b]; j < high; j += prime)
        sieve[j - low] = 0;
      next[b] = j;
    }

    int64_t count_low_high = 0;
    for (int64_t i = low; i < high; i++)
      count_low_high += sieve[i - low];

    // For k + 1 <= b <= pi_x_star
    // Find all special leaves: n = primes[b] * m
    // In the interval: low <= (x / n) < high
    // Which satisfy:  mu[m] != 0 && lpf[m] > primes[b] && mpf[m] <= y
    for (; b <= pi_x_star; b++)
    {
      int64_t prime = primes[b];
      int64_t max_m = min(x / (prime * low), x / ipow(prime, 3), z);
      int64_t min_m = max(x / (prime * high), z / prime, prime);
      int64_t count = 0;
      int64_t i = 0;

      if (prime >= max_m)
        break;

      for (int64_t m = max_m; m > min_m; m--)
      {
        if (mu[m] != 0 &&
            lpf[m] > prime &&
            mpf[m] <= y)
        {
          // We have found a special leaf. Compute it's contribution
          // phi(x / (primes[b] * m), b - 1) by counting the number
          // of unsieved elements <= x / (primes[b] * m) after having
          // removed the multiples of the first b - 1 primes.
          int64_t xpm = x / (prime * m);
          int64_t stop = xpm - low;
          for (; i <= stop; i++)
            count += sieve[i];
          int64_t phi_xpm = phi[b] + count;
          sum -= mu[m] * phi_xpm;
        }
      }

      phi[b] += count_low_high;

      // Remove the multiples of the b-th prime
      int64_t j = next[b];
      for (; j < high; j += prime * 2)
      {
        count_low_high -= sieve[j - low];
        sieve[j - low] = 0;
      }
      next[b] = j;
    }
  }

  print("D", sum, time);
  return sum;
}

} // namespace
