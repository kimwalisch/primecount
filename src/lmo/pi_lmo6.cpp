///
/// @file  pi_lmo6.cpp
/// @brief Implementation of the Lagarias-Miller-Odlyzko prime counting
///        algorithm. This version uses the modified algorithm as
///        described in section 5 (pages 556-557) in the paper
///        "Computing pi(x) The Meissel-Lehmer Method", Mathematics of
///        Computation, 44 (1985), by J. C. Lagarias, V. S. Miller and
///        A. M. Odlyzko.
///
///        In this version the special leaves for c + 1 <= b < pi_sqrty
///        have been split up into 2 categories:
///           1) The special leaves that are a product of 2 primes.
///           2) The special leaves that are a product of a prime and a
///              square free integer (which must not be prime).
///        Although this split up can give up to 15 percent speed
///        improvement it uses considerably more memory.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <pmath.hpp>
#include <pi_bsearch.hpp>
#include <PhiTiny.hpp>
#include <tos_counters.hpp>

#include <stdint.h>
#include <algorithm>
#include <vector>

using namespace std;

namespace {

/// Cross-off the multiples of prime in the sieve array.
/// For each element that is unmarked the first time update
/// the special counters tree data structure.
///
template <typename T1, typename T2>
void cross_off(int64_t prime, int64_t low, int64_t high, int64_t& next_multiple, T1& sieve, T2& counters)
{
  int64_t segment_size = sieve.size();
  int64_t k = next_multiple;

  for (; k < high; k += prime * 2)
  {
    if (sieve[k - low])
    {
      sieve[k - low] = 0;
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
           int64_t pi_y,
           int64_t c,
           vector<int32_t>& primes,
           vector<int32_t>& lpf,
           vector<int32_t>& mu)
{
  int64_t limit = x / y + 1;
  int64_t segment_size = next_power_of_2(isqrt(limit));
  int64_t sqrty = isqrt(y);
  int64_t pi_sqrty = pi_bsearch(primes, sqrty);
  int64_t S2_result = 0;

  vector<char> sieve(segment_size);
  vector<int32_t> counters(segment_size);
  vector<int32_t> kim(pi_sqrty, 0);
  vector<int32_t> pi = make_pi(y);
  vector<int64_t> next(primes.begin(), primes.end());
  vector<int64_t> phi(primes.size(), 0);

  vector<vector<int32_t> > square_free_candidates = 
      generate_square_free_candidates(c, y, lpf, mu, pi, primes);

  vector<vector<int32_t>::reverse_iterator > square_free_iters(pi_sqrty);
  for (int64_t i = 0; i < pi_sqrty; i++)
    square_free_iters[i] = square_free_candidates[i].rbegin();

  // Segmented sieve of Eratosthenes
  for (int64_t low = 1; low < limit; low += segment_size)
  {
    fill(sieve.begin(), sieve.end(), 1);

    // Current segment = interval [low, high[
    int64_t high = min(low + segment_size, limit);
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

    // For c + 1 <= b < pi_sqrty
    // Find all special leaves: n = primes[b] * m
    // which satisfy:  mu[m] != 0 && primes[b] < lpf[m], low <= (x / n) < high
    for (; b < pi_sqrty; b++)
    {
      int64_t prime = primes[b];
      int64_t min_m = max(x / (prime * high), y / prime);
      int64_t l = pi[min(x / (prime * low), y)];
      min_m = in_between(prime, min_m, y);
      int64_t l_min = pi[min_m];
      vector<int32_t>::reverse_iterator iter = square_free_iters[b];

      if (prime >= primes[l])
        goto next_segment;

      // Special leaves which are a product of 2 primes
      for (; l > l_min; l--)
      {
        int64_t n = prime * primes[l];
        int64_t count = cnt_query(counters, (x / n) - low);
        int64_t phi_xn = phi[b] + count;
        S2_result += phi_xn;
      }

      // Special leaves which are a product of a prime and a
      // square free integer which must satisfy:
      // !is_prime(square_free) && prime < least_prime_factor[square_free]
      for (; *iter > min_m; iter++)
      {
        int64_t square_free = *iter;
        int64_t n = prime * square_free;
        int64_t count = cnt_query(counters, (x / n) - low);
        int64_t phi_xn = phi[b] + count;
        S2_result -= mu[square_free] * phi_xn;
      }

      square_free_iters[b] = iter;
      phi[b] += cnt_query(counters, (high - 1) - low);
      cross_off(prime, low, high, next[b], sieve, counters);
    }

    // For pi_sqrty <= b < pi_y
    // Find all special leaves: n = primes[b] * prime2
    // which satisfy: low <= (x / n) < high
    for (; b < pi_y; b++)
    {
      int64_t prime = primes[b];
      int64_t l = pi[min(x / (prime * low), y)];
      int64_t min_m = max(x / (prime * high), y / prime);
      min_m = in_between(prime, min_m, y);
      int64_t min_l = pi[min_m];

      if (prime >= primes[l])
        goto next_segment;

      for (; l > min_l; l--)
      {
        int64_t n = prime * primes[l];
        int64_t count = cnt_query(counters, (x / n) - low);
        int64_t phi_xn = phi[b] + count;
        S2_result += phi_xn;
      }

      phi[b] += cnt_query(counters, (high - 1) - low);
      cross_off(prime, low, high, next[b], sieve, counters);
    }

    next_segment:;
  }

  return S2_result;
}

} // namespace

namespace primecount {

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo6(int64_t x)
{
  if (x < 2)
    return 0;

  // Optimization factor, see:
  // J. C. Lagarias, V. S. Miller, and A. M. Odlyzko, Computing pi(x): The Meissel-
  // Lehmer method, Mathematics of Computation, 44 (1985), p. 556.
  double beta = 1.0;
  double alpha = in_between(1, log(log((double) x)) * beta, iroot<6>(x));
  int64_t x13 = iroot<3>(x);
  int64_t y = (int64_t)(x13 * alpha);

  vector<int32_t> mu = make_moebius(y);
  vector<int32_t> lpf = make_least_prime_factor(y);
  vector<int32_t> primes;
  primes.push_back(0);
  primesieve::generate_primes(y, &primes);

  int64_t pi_y = primes.size() - 1;
  int64_t c = min<int64_t>(PhiTiny::MAX_A, pi_y);
  int64_t s1 = S1(x, y, c, primes, lpf , mu);
  int64_t s2 = S2(x, y, pi_y, c, primes, lpf , mu);
  int64_t p2 = P2(x, y, 1);
  int64_t phi = s1 + s2;
  int64_t sum = phi + pi_y - 1 - p2;

  return sum;
}

} // namespace primecount
