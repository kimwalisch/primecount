///
/// @file  PhiCache.cpp
/// @brief The PhiCache class calculates phi(x, a) using the recursive
///        formula: phi(x, a) = phi(x, a - 1) - phi(x / primes[a], a - 1).
///        I have added a cache to my implementation in which results
///        of phi(x, a) are stored if x < 2^16 and a <= 500. The cache
///        speeds up the calculations by at least 3 orders of
///        magnitude near 10^15.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PhiCache.hpp>
#include <PhiTiny.hpp>
#include <pmath.hpp>

#include <stdint.h>
#include <algorithm>
#include <limits>
#include <vector>
#include <cassert>

using namespace std;

namespace {

inline int64_t fast_div(int64_t x, int32_t y)
{
  // Avoid slow 64-bit division
  if (x <= numeric_limits<uint32_t>::max())
    return ((uint32_t) x) / y;
  return x / y;
}

} // namespace

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

  if (x < primes_[a])
    sum = sign;
  else if (is_phi_tiny(a))
    sum = phi_tiny(x, a) * sign;
  else if (is_phi_bsearch(x, a))
    sum = phi_bsearch(x, a) * sign;
  else
  {
    int64_t iters = pi_bsearch(primes_, a, isqrt(x));
    int64_t c = min(iters, PhiTiny::max_a());
    sum = (a - iters) * -sign;
    sum += phi_tiny(x, c) * sign;

    for (int64_t a2 = c; a2 < iters; a2++)
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
