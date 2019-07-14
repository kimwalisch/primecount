///
/// @file  pi_deleglise_rivat1.cpp
/// @brief Simple demonstration implementation of the Deleglise-Rivat
///        prime counting algorithm. The Deleglise-Rivat algorithm is
///        an improvement over the Lagarias-Miller-Odlyzko algorithm,
///        in the Deleglise-Rivat algorithm the special leaves
///        S2(x, a) have been split up into trivial special leaves,
///        easy special leaves and hard special leaves.
///
///        Deleglise-Rivat formula:
///        pi(x) = pi(y) + S1(x, a) + S2(x, a) - 1 - P2(x, a)
///        S2(x, a) = S2_trivial(x, a) + S2_easy(x, a) + S2_hard(x, a)
///        with y = alpha * x^(1/3), a = pi(y)
///
///        This implementation is based on the paper:
///        Tomás Oliveira e Silva, Computing pi(x): the combinatorial
///        method, Revista do DETUA, vol. 4, no. 6, March 2006,
///        pp. 759-768.
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <BinaryIndexedTree.hpp>
#include <generate.hpp>
#include <min.hpp>
#include <imath.hpp>
#include <PhiTiny.hpp>
#include <S1.hpp>
#include <S2.hpp>

#include <stdint.h>
#include <algorithm>
#include <vector>

using namespace std;
using namespace primecount;

namespace {

/// Cross-off the multiples of prime in the sieve array.
/// For each element that is unmarked the first time update
/// the binary indexed tree data structure.
///
template <typename Sieve>
void cross_off(int64_t prime,
               int64_t low,
               int64_t high,
               int64_t& multiple,
               Sieve& sieve,
               BinaryIndexedTree& tree)
{
  int64_t m = multiple;

  for (; m < high; m += prime * 2)
  {
    if (sieve[m - low])
    {
      sieve[m - low] = 0;
      tree.update(m - low);
    }
  }

  multiple = m;
}

/// Compute the contribution of the hard special leaves using a
/// segmented sieve in order to reduce memory usage.
///
/// Note that in the Deleglise-Rivat paper it is suggested to use a
/// segment size of y. In practice however this uses too much memory
/// especially when using multi-threading. Hence we are using a
/// segment size of sqrt(z) as suggested in Xavier Gourdon's paper.
/// In primecount's implementation a segment size of sqrt(z) seems
/// ideal since slightly increasing the segment size decreases
/// performance because of cache misses and slightly decreasing the
/// segment size also decreases performance.
///
int64_t S2_hard(int64_t x,
                int64_t y,
                int64_t z,
                int64_t c,
                vector<int32_t>& pi,
                vector<int32_t>& lpf,
                vector<int32_t>& mu)
{
  int64_t limit = z + 1;
  int64_t segment_size = next_power_of_2(isqrt(limit));
  int64_t pi_sqrty = pi[isqrt(y)];
  int64_t pi_sqrtz = pi[min(isqrt(z), y)];
  int64_t s2_hard = 0;

  auto primes = generate_primes<int32_t>(y);

  vector<char> sieve(segment_size);
  vector<int64_t> next(primes.begin(), primes.end());
  vector<int64_t> phi(primes.size(), 0);
  BinaryIndexedTree tree;

  // segmented sieve of Eratosthenes
  for (int64_t low = 1; low < limit; low += segment_size)
  {
    // current segment [low, high[
    int64_t high = min(low + segment_size, limit);
    int64_t b = 1;

    fill(sieve.begin(), sieve.end(), 1);

    // pre-sieve multiples of first c primes
    for (; b <= c; b++)
    {
      int64_t k = next[b];
      for (int64_t prime = primes[b]; k < high; k += prime)
        sieve[k - low] = 0;
      next[b] = k;
    }

    // initialize binary indexed tree from sieve
    tree.init(sieve);

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

      for (int64_t m = max_m; m > min_m; m--)
      {
        if (mu[m] != 0 && prime < lpf[m])
        {
          int64_t n = prime * m;
          int64_t count = tree.count(low, x / n);
          int64_t phi_xn = phi[b] + count;
          s2_hard -= mu[m] * phi_xn;
        }
      }

      phi[b] += tree.count(low, high - 1);
      cross_off(prime, low, high, next[b], sieve, tree);
    }

    // For pi_sqrty < b <= pi_sqrtz
    // Find all hard special leaves: n = primes[b] * primes[l]
    // which satisfy: low <= (x / n) < high
    for (; b <= pi_sqrtz; b++)
    {
      int64_t prime = primes[b];
      int64_t l = pi[min(x / (prime * low), z / prime, y)];
      int64_t min_hard = max(x / (prime * high), prime);

      if (prime >= primes[l])
        goto next_segment;

      for (; primes[l] > min_hard; l--)
      {
        int64_t n = prime * primes[l];
        int64_t xn = x / n;
        int64_t count = tree.count(low, xn);
        int64_t phi_xn = phi[b] + count;
        s2_hard += phi_xn;
      }

      phi[b] += tree.count(low, high - 1);
      cross_off(prime, low, high, next[b], sieve, tree);
    }

    next_segment:;
  }

  return s2_hard;
}

/// Calculate the contribution of the special leaves
int64_t S2(int64_t x,
           int64_t y,
           int64_t z,
           int64_t c,
           vector<int32_t>& lpf,
           vector<int32_t>& mu)
{
  auto pi = generate_pi(y);

  int64_t s2_trivial = S2_trivial(x, y, z, c);
  int64_t s2_easy = S2_easy(x, y, z, c, 1);
  int64_t s2_hard = S2_hard(x, y, z, c, pi, lpf, mu);
  int64_t s2 = s2_trivial + s2_easy + s2_hard;

  return s2;
}

} // namespace

namespace primecount {

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2)
/// Memory usage: O(x^(1/3) * (log x)^3)
///
int64_t pi_deleglise_rivat1(int64_t x)
{
  if (x < 2)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  int64_t x13 = iroot<3>(x);
  int64_t y = (int64_t) (x13 * alpha);
  int64_t z = x / y;
  int64_t c = PhiTiny::get_c(y);
  int64_t p2 = P2(x, y, 1);

  auto mu = generate_moebius(y);
  auto lpf = generate_lpf(y);

  int64_t pi_y = pi_legendre(y);
  int64_t s1 = S1(x, y, c, 1);
  int64_t s2 = S2(x, y, z, c, lpf, mu);
  int64_t phi = s1 + s2;
  int64_t sum = phi + pi_y - 1 - p2;

  return sum;
}

} // namespace
