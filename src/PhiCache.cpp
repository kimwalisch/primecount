///
/// @file  PhiCache.cpp
/// @brief The PhiCache class calculates the partial sieve function
///        (a.k.a. Legendre-sum) using the recursive formula:
///        phi(x, a) = phi(x, a - 1) - phi(x / primes[a], a - 1).
///        phi(x, a) counts the numbers <= x that are not divisible by
///        any of the first a primes. The algorithm used is an
///        optimized version of the algorithm described in Tomás
///        Oliveira e Silva's paper [1]. I have added 5 optimizations
///        to my implementation which significantly speed up the
///        calculation:
///
///        * Cache results of phi(x, a)
///        * Calculate phi(x, a) using formula [2] if a <= 6
///        * Calculate phi(x, a) using binary search
///        * Calculate all phi(x, a) = 1 upfront
///        * Stop recursion at c instead of 1
///
///       [1] Tomás Oliveira e Silva, Computing pi(x): the combinatorial
///           method, Revista do DETUA, vol. 4, no. 6, March 2006, p. 761.
///           http://sweet.ua.pt/tos/bib/5.4.pdf
///       [2] phi(x, a) = (x / pp) * φ(pp) + phi(x % pp, a)
///           with pp = 2 * 3 * ... * prime[a] 
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PhiCache.hpp>
#include <PhiTiny.hpp>
#include <pmath.hpp>
#include <fast_div.hpp>

#include <stdint.h>
#include <algorithm>
#include <limits>
#include <vector>
#include <cassert>

using namespace std;

namespace primecount {

PhiCache::PhiCache(const vector<int32_t>& primes) :
  primes_(primes),
  bytes_(0)
{
  // primecount uses 1-indexing i.e. primes[1] = 2
  assert(primes_[0] == 0);
  size_t max_size = CACHE_A_LIMIT + 1;
  cache_.resize(min(primes.size(), max_size));
}

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
///
int64_t PhiCache::phi(int64_t x, int64_t a)
{
  if (x < 1) return 0;
  if (a > x) return 1;
  if (a < 1) return x;

  if (primes_.at(a) >= x)
    return 1;

  return phi(x, a, 1);
}

/// Calculate phi(x, a) using the recursive formula:
/// phi(x, a) = phi(x, a - 1) - phi(x / primes_[a], a - 1)
///
int64_t PhiCache::phi(int64_t x, int64_t a, int sign)
{
  int64_t sum;

  if (x <= primes_[a])
    sum = sign;
  else if (is_phi_tiny(a))
    sum = phi_tiny(x, a) * sign;
  else if (is_phi_bsearch(x, a))
    sum = phi_bsearch(x, a) * sign;
  else
  {
    // Move out of the loop the calculations where phi(x2, a2) = 1
    // phi(x, a) = 1 if primes_[a] >= x
    // x2 = x / primes_[a2 + 1]
    // phi(x2, a2) = 1 if primes_[a2] >= x / primes_[a2 + 1]
    // phi(x2, a2) = 1 if primes_[a2] >= sqrt(x)
    // phi(x2, a2) = 1 if a2 >= pi(sqrt(x))
    // \sum_{a2 = pi(sqrt(x))}^{a - 1} phi(x2, a2) = a - pi(sqrt(x))
    //
    int64_t pi_sqrtx = pi_bsearch(primes_, a, isqrt(x));
    sum = (a - pi_sqrtx) * -sign;

    // phi(x, 1) - \sum_{a2 = 1}^{c - 1} phi(x / primes_[a2 + 1], a2) = phi(x, c)
    int64_t c = min(pi_sqrtx, PhiTiny::max_a());
    sum += phi_tiny(x, c) * sign;

    for (int64_t a2 = c; a2 < pi_sqrtx; a2++)
    {
      int64_t x2 = fast_div(x, primes_[a2 + 1]);
      if (is_cached(x2, a2))
        sum += cache_[a2][x2] * -sign;
      else
        sum += phi(x2, a2, -sign);
    }
  }

  if (write_to_cache(x, a))
    cache_[a][x] = (uint16_t) (sum * sign);

  return sum;
}

/// Binary search phi(x, a)
int64_t PhiCache::phi_bsearch(int64_t x, int64_t a) const
{
  int64_t pix = pi_bsearch(primes_, x);
  return pix - a + 1;
}

bool PhiCache::is_phi_bsearch(int64_t x, int64_t a) const
{
  return x <= primes_.back() &&
         x < isquare(primes_[a + 1]);
}

bool PhiCache::write_to_cache(int64_t x, int64_t a)
{
  if (a > CACHE_A_LIMIT || x > numeric_limits<uint16_t>::max())
    return false;

  if (x >= cache_size(a))
  {
    if (bytes_ > CACHE_BYTES_LIMIT)
      return false;
    bytes_ += (x + 1 - cache_size(a)) * 2;
    cache_[a].resize(x + 1, 0);
  }

  return true;
}

} // namespace primecount
