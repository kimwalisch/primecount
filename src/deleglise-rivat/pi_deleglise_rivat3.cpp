///
/// @file  pi_deleglise_rivat3.cpp
/// @brief Implementation of the Lagarias-Miller-Odlyzko prime counting
///        algorithm with the improvements of Deleglise and Rivat.
///        This version uses compression (see FactorTable& PiTable) to
///        reduce the memory usage.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <FactorTable.hpp>
#include <primecount-internal.hpp>
#include <BitSieve.hpp>
#include <generate.hpp>
#include <pmath.hpp>
#include <PhiTiny.hpp>
#include <S1.hpp>
#include <tos_counters.hpp>

#include <stdint.h>
#include <algorithm>
#include <iostream>
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
  int64_t k = next_multiple;

  for (; k < high; k += prime * 2)
  {
    if (sieve[k - low])
    {
      sieve.unset(k - low);
      cnt_update(counters, k - low, segment_size);
    }
  }
  next_multiple = k;
}

/// Calculate the contribution of the special leaves.
/// @see ../docs/computing-special-leaves.md
/// @pre y > 0 && c > 1
///
int64_t S2(int64_t x,
           int64_t y,
           int64_t z,
           int64_t c,
           vector<int32_t>& primes,
           FactorTable<uint16_t>& factors)
{
  PiTable pi(y);
  int64_t pi_y = pi(y);
  int64_t pi_sqrty = pi(isqrt(y));
  int64_t pi_sqrtz = pi(min(isqrt(z), y));
  int64_t limit = z + 1;
  int64_t segment_size = next_power_of_2(isqrt(limit));
  int64_t S2_result = 0;

  BitSieve sieve(segment_size);
  vector<int32_t> counters(segment_size);
  vector<int64_t> next(primes.begin(), primes.begin() + pi_sqrtz + 1);
  vector<int64_t> phi(pi_sqrtz + 1, 0);

  // Segmented sieve of Eratosthenes
  for (int64_t low = 1; low < limit; low += segment_size)
  {
    // Current segment = interval [low, high[
    int64_t high = min(low + segment_size, limit);
    int64_t b = c + 1;

    // check if we need the sieve
    if (c < pi_sqrtz)
    {
      sieve.memset(low);

      // phi(y, i) nodes with i <= c do not contribute to S2, so we
      // simply sieve out the multiples of the first c primes
      for (int64_t i = 2; i <= c; i++)
      {
        int64_t k = next[i];
        for (int64_t prime = primes[i]; k < high; k += prime * 2)
          sieve.unset(k - low);
        next[i] = k;
      }

      // Initialize special tree data structure from sieve
      cnt_finit(sieve, counters, segment_size);
    }

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

    // For pi_sqrty <= b < pi_y
    // Find all special leaves: n = primes[b] * primes[l]
    // which satisfy: low <= (x / n) < high
    for (; b < pi_y; b++)
    {
      int64_t prime = primes[b];
      int64_t l = pi(min(x / (prime * low), y));

      if (prime >= primes[l])
        goto next_segment;

      int64_t min_hard_leaf = max3(x / (prime * high), y / prime, prime);
      int64_t min_trivial_leaf = min(x / (prime * prime), y);
      int64_t min_clustered_easy_leaf = min(isqrt(x / prime), y);
      int64_t min_sparse_easy_leaf = min(z / prime, y);

      min_trivial_leaf = max(min_hard_leaf, min_trivial_leaf);
      min_clustered_easy_leaf = max(min_hard_leaf, min_clustered_easy_leaf);
      min_sparse_easy_leaf = max(min_hard_leaf, min_sparse_easy_leaf);

      // Find all trivial leaves which satisfy:
      // phi(x / (primes[b] * primes[l]), b - 1) = 1
      if (primes[l] > min_trivial_leaf)
      {
        int64_t l_min = pi(min_trivial_leaf);
        S2_result += l - l_min;
        l = l_min;
      }

      // Find all clustered easy leaves which satisfy:
      // x / n <= y such that phi(x / n, b - 1) = pi(x / n) - b + 2
      // And phi(x / n, b - 1) == phi(x / m, b - 1)
      while (primes[l] > min_clustered_easy_leaf)
      {
        int64_t n = prime * primes[l];
        int64_t xn = x / n;
        assert(xn < isquare(primes[b]));
        int64_t phi_xn = pi(xn) - b + 2;
        int64_t m = prime * primes[b + phi_xn - 1];
        int64_t xm = max(x / m, min_clustered_easy_leaf);
        int64_t l2 = pi(xm);
        S2_result += phi_xn * (l - l2);
        l = l2;
      }

      // Find all sparse easy leaves which satisfy:
      // x / n <= y such that phi(x / n, b - 1) = pi(x / n) - b + 2
      for (; primes[l] > min_sparse_easy_leaf; l--)
      {
        int64_t n = prime * primes[l];
        int64_t xn = x / n;
        assert(xn < isquare(primes[b]));
        S2_result += pi(xn) - b + 2;
      }

      if (b <= pi_sqrtz)
      {
        // Find all hard leaves which satisfy:
        // low <= (x / n) < high
        for (; primes[l] > min_hard_leaf; l--)
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
    }

    next_segment:;
  }

  return S2_result;
}

/// alpha is a tuning factor which should grow like (log(x))^3
/// for the Deleglise-Rivat prime counting algorithm.
///
double compute_alpha(int64_t x)
{
  double d = (double) x;
  double alpha = log(d) * log(d) * log(d) / 1200;
  return in_between(1, alpha, iroot<6>(x));
}

} // namespace

namespace primecount {

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
int64_t pi_deleglise_rivat3(int64_t x)
{
  if (x < 2)
    return 0;

  double alpha = compute_alpha(x);
  int64_t y = (int64_t) (alpha * iroot<3>(x));
  int64_t z = x / y;

  if (print_status())
  {
    cout << endl;
    cout << "=== pi_deleglise_rivat3 ===" << endl;
    cout << "pi(x) = S1 + S2 - 1 - P2" << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
    cout << "threads = " << "1" << endl;
  }

  int64_t p2 = P2(x, y, 1);

  vector<int32_t> primes = generate_primes(y);
  FactorTable<uint16_t> factors(y);

  int64_t pi_y = pi_bsearch(primes, y);
  int64_t c = min(pi_y, PhiTiny::max_a());
  int64_t s1 = S1(x, y, c, primes[c], factors, 1);
  int64_t s2 = S2(x, y, z, c, primes, factors);
  int64_t phi = s1 + s2;
  int64_t sum = phi + pi_y - 1 - p2;

  if (print_status())
    cout << endl;

  return sum;
}

} // namespace primecount
