///
/// @file  pi_lmo4.cpp
/// @brief Simple demonstration implementation of the
///        Lagarias-Miller-Odlyzko prime counting algorithm.
///        This implementation uses the segmented sieve of
///        Eratosthenes to calculate the contribution of the special
///        leaves and a special tree data structure (a.k.a. Fenwick
///        tree) to count the number of unsieved elements.
///
///        Lagarias-Miller-Odlyzko formula:
///        pi(x) = pi(y) + S1(x, a) + S2(x, a) - 1 - P2(x, a)
///        with y = x^(1/3), a = pi(y)
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "BinaryIndexedTree.hpp"

#include <primecount-internal.hpp>
#include <generate_primes.hpp>
#include <imath.hpp>
#include <min.hpp>
#include <PhiTiny.hpp>
#include <Vector.hpp>
#include <S.hpp>

#include <stdint.h>
#include <algorithm>

using namespace primecount;

namespace {

/// Cross-off the multiples of prime in the sieve array.
/// For each element that is unmarked the first time update
/// the binary indexed tree data structure.
///
template <typename S, typename T>
void cross_off(int64_t prime,
               int64_t low,
               int64_t high,
               int64_t& next_multiple,
               S& sieve,
               T& tree)
{
  int64_t m = next_multiple;

  for (; m < high; m += prime * 2)
  {
    if (sieve[m - low])
    {
      sieve[m - low] = 0;
      tree.update(m - low);
    }
  }

  next_multiple = m;
}

/// Calculate the contribution of the special leaves
int64_t S2(int64_t x,
           int64_t y,
           int64_t c,
           int64_t pi_y,
           const Vector<int32_t>& primes,
           const Vector<int32_t>& lpf,
           const Vector<int32_t>& mu)
{
  int64_t limit = x / y;
  int64_t segment_size = next_power_of_2(isqrt(limit));
  int64_t s2 = 0;

  BinaryIndexedTree tree;
  Vector<bool> sieve(segment_size);
  Vector<int64_t> next(primes.size());
  Vector<int64_t> phi(primes.size());
  std::copy(primes.begin(), primes.end(), next.begin());
  std::fill(phi.begin(), phi.end(), 0);

  // segmented sieve of Eratosthenes
  for (int64_t low = 1; low < limit; low += segment_size)
  {
    // current segment [low, high[
    int64_t high = min(low + segment_size, limit);

    std::fill(sieve.begin(), sieve.end(), 1);

    // phi(y, b) nodes with b <= c do not contribute to S2, so
    // we sieve out the multiples of the first c primes
    for (int64_t b = 1; b <= c; b++)
    {
      int64_t k = next[b];
      for (int64_t prime = primes[b]; k < high; k += prime)
        sieve[k - low] = 0;
      next[b] = k;
    }

    // initialize binary indexed tree from sieve
    tree.init(sieve);

    for (int64_t b = c + 1; b < pi_y; b++)
    {
      int64_t prime = primes[b];
      int64_t min_m = max(x / (prime * high), y / prime);
      int64_t max_m = min(x / (prime * low), y);

      // Obviously if (prime >= max_m) then (prime >= lpf[max_m])
      // hence (prime < lpf[m]) will always evaluate to
      // false and no special leaves are possible.
      if (prime >= max_m)
        break;

      for (int64_t m = max_m; m > min_m; m--)
      {
        if (mu[m] != 0 && prime < lpf[m])
        {
          // We have found a special leaf. Compute it's contribution
          // phi(x / (primes[b] * m), b - 1) by counting the number
          // of unsieved elements <= x / (primes[b] * m) after having
          // removed the multiples of the first b - 1 primes.
          //
          int64_t n = prime * m;
          int64_t count = tree.count(low, x / n);
          int64_t phi_xn = phi[b] + count;
          s2 -= mu[m] * phi_xn;
        }
      }

      // save the number of unsieved elements
      phi[b] += tree.count(low, high - 1);

      // remove the multiples of the b-th prime
      cross_off(prime, low, high, next[b], sieve, tree);
    }
  }

  return s2;
}

} // namespace

namespace primecount {

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3))
/// Memory usage: O(x^(1/3) * (log x)^2)
///
int64_t pi_lmo4(int64_t x)
{
  if (x < 2)
    return 0;

  int threads = 1;
  double alpha = get_alpha_lmo(x);
  int64_t x13 = iroot<3>(x);
  int64_t y = (int64_t) (x13 * alpha);
  int64_t c = PhiTiny::get_c(y);

  auto primes = generate_primes<int32_t>(y);
  auto lpf = generate_lpf(y);
  auto mu = generate_moebius(y);

  int64_t pi_y = primes.size() - 1;
  int64_t p2 = P2(x, y, pi_y, threads);
  int64_t s1 = S1(x, y, c, threads);
  int64_t s2 = S2(x, y, c, pi_y, primes, lpf, mu);
  int64_t phi = s1 + s2;
  int64_t sum = phi + pi_y - 1 - p2;

  return sum;
}

} // namespace
