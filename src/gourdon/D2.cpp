///
/// @file  D2.cpp
/// @brief Simple demonstration implementation of the D(x, y) formula
///        in Xavier Gourdon's prime counting algorithm. This
///        implementation runs single threaded and does not use the
///        highly optimized Sieve.cpp.
///
///        In this implementation the hard special leaves have been
///        split up into 2 distinct types. Below sqrt(z) the leaves
///        are composed of a prime and a square free number. But when
///        the prime factors are > sqrt(z) then all leaves are
///        composed of exactly 2 primes.
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
#include <PiTable.hpp>

#include <stdint.h>
#include <vector>

using namespace std;

namespace {

/// Remove the multiples of the b-th prime from the sieve array.
/// @return: Count of numbers unset for the 1st time.
///
template <typename Sieve>
int64_t cross_off(Sieve& sieve,
                  int64_t prime,
                  int64_t& next_multiple,
                  int64_t low,
                  int64_t high)
{
  int64_t n = next_multiple;
  int64_t count = 0;

  for (; n < high; n += prime * 2)
  {
    count += sieve[n - low];
    sieve[n - low] = 0;
  }

  next_multiple = n;
  return count;
}

} // namespace

namespace primecount {

int64_t D(int64_t x,
          int64_t y,
          int64_t z,
          int64_t k)
{
  print("");
  print("=== D(x, y) ===");
  print(x, y, z, k, 1);

  double time = get_time();
  int64_t sum = 0;
  int64_t limit = x / z + 1;
  int64_t segment_size = isqrt(limit);
  int64_t x_star = get_x_star_gourdon(x, y);

  PiTable pi(y);
  auto primes = generate_primes<int32_t>(y);

  auto mu = generate_moebius(z);
  auto lpf = generate_lpf(z);
  auto mpf = generate_mpf(z);

  int64_t pi_sqrtz = pi[isqrt(z)];
  int64_t pi_x_star = pi[x_star];

  vector<char> sieve(segment_size);
  vector<int64_t> phi(pi_x_star + 1, 0);
  vector<int64_t> next;

  next.reserve(pi_x_star + 1);
  for (int64_t i = 0; i <= pi_x_star; i++)
    next.push_back(primes[i]);

  // Segmented sieve of Eratosthenes
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

    // For k + 1 <= b <= pi_sqrtz
    // Find all special leaves: n = primes[b] * m
    // In the interval: low <= (x / n) < high
    // Which satisfy:  mu[m] != 0 && lpf[m] > primes[b] && mpf[m] <= y
    for (; b <= pi_sqrtz; b++)
    {
      int64_t prime = primes[b];
      int64_t max_m = min(x / (prime * low), x / ipow(prime, 3), z);
      int64_t min_m = max(x / (prime * high), z / prime, prime);
      int64_t count = 0;
      int64_t i = 0;

      if (prime >= max_m)
        goto next_segment;

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
      count_low_high -= cross_off(sieve, prime, next[b], low, high);
    }

    // For pi_sqrtz < b <= pi_x_star
    // Find all special leaves: n = primes[b] * prime2
    // which satisfy: low <= (x / n) < high && prime2 <= y
    for (; b <= pi_x_star; b++)
    {
      int64_t prime = primes[b];
      int64_t max_m = min(x / (prime * low), x / ipow(prime, 3), y);
      int64_t min_m = max(x / (prime * high), z / prime, prime);
      int64_t l = pi[max_m];
      int64_t i = 0;
      int64_t count = 0;

      if (prime >= primes[l])
        goto next_segment;

      for (; primes[l] > min_m; l--)
      {
        int64_t xpq = x / (prime * primes[l]);
        int64_t stop = xpq - low;
        for (; i <= stop; i++)
          count += sieve[i];
        int64_t phi_xpq = phi[b] + count;
        sum += phi_xpq;
      }

      phi[b] += count_low_high;
      count_low_high -= cross_off(sieve, prime, next[b], low, high);
    }

    next_segment:;
  }

  print("D", sum, time);
  return sum;
}

} // namespace
