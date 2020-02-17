///
/// @file  D3.cpp
/// @brief Single threaded implementation of the D(x, y) formula in
///        Xavier Gourdon's prime counting algorithm. This
///        implementation uses the highly optimized Sieve class.
///
///        This implementation also uses the DFactorTable lookup
///        table instead of the mu, lpf and mpf lookup tables.
///        DFactorTable uses much less memory and allows to check
///        more quickly whether a number is a leaf or not.
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
#include <Sieve.hpp>

#include "DFactorTable.hpp"

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
  print_gourdon_vars(x, y, z, k, 1);

  double time = get_time();
  int64_t sum = 0;
  int64_t limit = x / z;
  int64_t segment_size = Sieve::get_segment_size(isqrt(limit));
  int64_t x_star = get_x_star_gourdon(x, y);
  int64_t low = 0;

  auto primes = generate_primes<int32_t>(y);
  DFactorTable<uint16_t> factor(y, z, 1);
  Sieve sieve(low, segment_size, primes.size());

  PiTable pi(y);
  int64_t pi_sqrtz = pi[isqrt(z)];
  int64_t pi_x_star = pi[x_star];
  vector<int64_t> phi(pi_x_star + 1, 0);

  // Segmented sieve of Eratosthenes
  for (; low < limit; low += segment_size)
  {
    // Current segment [low, high[
    int64_t high = min(low + segment_size, limit);
    int64_t low1 = max(low, 1);

    sieve.pre_sieve(primes, k, low, high);
    int64_t b = k + 1;

    // For k + 1 <= b <= pi_sqrtz
    // Find all special leaves in the current segment that are
    // composed of a prime and a square free number:
    // low <= x / (primes[b] * m) < high
    for (; b <= pi_sqrtz; b++)
    {
      int64_t prime = primes[b];
      int64_t max_m = min3(x / (prime * low1), x / ipow(prime, 3), z);
      int64_t min_m = max3(x / (prime * high), z / prime, prime);

      if (prime >= max_m)
        goto next_segment;

      min_m = factor.to_index(min_m);
      max_m = factor.to_index(max_m);

      for (int64_t m = max_m; m > min_m; m--)
      {
        // mu[m] != 0 && 
        // lpf[m] > prime && 
        // mpf[m] <= y
        if (prime < factor.is_leaf(m))
        {
          int64_t xpm = x / (prime * factor.to_number(m));
          int64_t stop = xpm - low;
          int64_t phi_xpm = phi[b] + sieve.count(stop);
          sum -= factor.mu(m) * phi_xpm;
        }
      }

      phi[b] += sieve.get_total_count();
      sieve.cross_off_count(prime, b);
    }

    // For pi_sqrtz < b <= pi_x_star
    // Find all special leaves in the current segment
    // that are composed of 2 primes:
    // low <= x / (primes[b] * primes[l]) < high
    for (; b <= pi_x_star; b++)
    {
      int64_t prime = primes[b];
      int64_t max_m = min3(x / (prime * low1), x / ipow(prime, 3), y);
      int64_t min_m = max3(x / (prime * high), z / prime, prime);
      int64_t l = pi[max_m];

      if (prime >= primes[l])
        goto next_segment;

      for (; primes[l] > min_m; l--)
      {
        int64_t xpq = x / (prime * primes[l]);
        int64_t stop = xpq - low;
        int64_t phi_xpq = phi[b] + sieve.count(stop);
        sum += phi_xpq;
      }

      phi[b] += sieve.get_total_count();
      sieve.cross_off_count(prime, b);
    }

    next_segment:;
  }

  print("D", sum, time);
  return sum;
}

} // namespace
