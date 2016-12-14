///
/// @file  pi_deleglise_rivat2.cpp
/// @brief Implementation of the Deleglise-Rivat prime counting
///        algorithm. Compared to pi_deleglise_rivat1.cpp this version
///        uses compression (FactorTable & PiTable) to reduce the
///        memory usage.
///
///        This implementation is based on the paper:
///        Tomás Oliveira e Silva, Computing pi(x): the combinatorial
///        method, Revista do DETUA, vol. 4, no. 6, March 2006,
///        pp. 759-768.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <FactorTable.hpp>
#include <primecount-internal.hpp>
#include <BitSieve.hpp>
#include <generate.hpp>
#include <min_max.hpp>
#include <imath.hpp>
#include <PhiTiny.hpp>
#include <tos_counters.hpp>
#include <S1.hpp>
#include "S2.hpp"

#include <stdint.h>
#include <algorithm>
#include <vector>

using namespace std;
using namespace primecount;

namespace {

/// Cross-off the multiples of prime in the sieve array.
/// For each element that is unmarked the first time update
/// the special counters tree data structure.
///
template <typename T>
void cross_off(int64_t prime,
               int64_t low,
               int64_t high,
               int64_t& next_multiple,
               BitSieve& sieve,
               T& counters)
{
  int64_t segment_size = sieve.size();
  int64_t m = next_multiple;

  for (; m < high; m += prime * 2)
  {
    if (sieve[m - low])
    {
      sieve.unset(m - low);
      cnt_update(counters, m - low, segment_size);
    }
  }

  next_multiple = m;
}

/// Calculate the contribution of the hard special leaves which
/// require use of a sieve (to reduce the memory usage).
///
int64_t S2_hard(int64_t x,
                int64_t y,
                int64_t z,
                int64_t c)
{
  print("");
  print("=== S2_hard(x, y) ===");
  print("Computation of the hard special leaves");

  PiTable pi(y);
  FactorTable<uint16_t> factors(y, 1);
  vector<int32_t> primes = generate_primes(y);

  int64_t limit = z + 1;
  int64_t segment_size = next_power_of_2(isqrt(limit));
  int64_t pi_sqrty = pi[isqrt(y)];
  int64_t pi_sqrtz = pi[min(isqrt(z), y)];
  int64_t S2_result = 0;
  double time = get_wtime();

  BitSieve sieve(segment_size);
  vector<int32_t> counters(segment_size);
  vector<int64_t> next(primes.begin(), primes.end());
  vector<int64_t> phi(primes.size(), 0);

  // Segmented sieve of Eratosthenes
  for (int64_t low = 1; low < limit; low += segment_size)
  {
    // Current segment = interval [low, high[
    int64_t high = min(low + segment_size, limit);
    int64_t b = c + 1;

    // pre-sieve the multiples of the first c primes
    sieve.pre_sieve(c, low);

    // Initialize special tree data structure from sieve
    cnt_finit(sieve, counters, segment_size);

    // For c + 1 <= b <= pi_sqrty
    // Find all special leaves: n = primes[b] * m, with mu[m] != 0 and primes[b] < lpf[m]
    // which satisfy: low <= (x / n) < high
    for (; b <= pi_sqrty; b++)
    {
      int64_t prime = primes[b];
      int64_t min_m = max(x / (prime * high), y / prime);
      int64_t max_m = min(x / (prime * low), y);

      if (prime >= max_m)
        goto next_segment;

      factors.to_index(&min_m);
      factors.to_index(&max_m);

      for (int64_t m = max_m; m > min_m; m--)
      {
        if (prime < factors.lpf(m))
        {
          int64_t n = prime * factors.get_number(m);
          int64_t count = cnt_query(counters, (x / n) - low);
          int64_t phi_xn = phi[b] + count;
          S2_result -= factors.mu(m) * phi_xn;
        }
      }

      phi[b] += cnt_query(counters, (high - 1) - low);
      cross_off(prime, low, high, next[b], sieve, counters);
    }

    // For pi_sqrty <= b <= pi_sqrtz
    // Find all hard special leaves: n = primes[b] * primes[l]
    // which satisfy: low <= (x / n) < high
    for (; b <= pi_sqrtz; b++)
    {
      int64_t prime = primes[b];
      int64_t l = pi[min3(x / (prime * low), z / prime, y)];
      int64_t min_hard = max3(x / (prime * high), y / prime, prime);

      if (prime >= primes[l])
        goto next_segment;

      for (; primes[l] > min_hard; l--)
      {
        int64_t n = prime * primes[l];
        int64_t xn = x / n;
        int64_t count = cnt_query(counters, xn - low);
        int64_t phi_xn = phi[b] + count;
        S2_result += phi_xn;
      }

      phi[b] += cnt_query(counters, (high - 1) - low);
      cross_off(prime, low, high, next[b], sieve, counters);
    }

    next_segment:;
  }

  print("S2_hard", S2_result, time);
  return S2_result;
}

/// Calculate the contribution of the special leaves.
/// @pre y > 0 && c > 1
///
int64_t S2(int64_t x,
           int64_t y,
           int64_t z,
           int64_t c)
{
  int64_t s2_trivial = S2_trivial(x, y, z, c, 1);
  int64_t s2_easy = S2_easy(x, y, z, c, 1);
  int64_t s2_hard = S2_hard(x, y, z, c);
  int64_t s2 = s2_trivial + s2_easy + s2_hard;

  return s2;
}

} // namespace

namespace primecount {

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
int64_t pi_deleglise_rivat2(int64_t x)
{
  if (x < 2)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  int64_t x13 = iroot<3>(x);
  int64_t y = (int64_t) (x13 * alpha);
  int64_t c = PhiTiny::get_c(y);
  int64_t z = x / y;

  print("");
  print("=== pi_deleglise_rivat2(x) ===");
  print("pi(x) = S1 + S2 + pi(y) - 1 - P2");
  print(x, y, z, c, alpha, 1);

  int64_t p2 = P2(x, y, 1);
  int64_t pi_y = pi_legendre(y, 1);
  int64_t s1 = S1(x, y, c, 1);
  int64_t s2 = S2(x, y, z, c);
  int64_t phi = s1 + s2;
  int64_t sum = phi + pi_y - 1 - p2;

  return sum;
}

} // namespace
