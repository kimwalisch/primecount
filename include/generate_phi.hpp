///
/// @file  generate_phi.hpp
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
///        * Calculate phi(x, a) in O(1) using formula [2] if a <= 6.
///        * Calculate phi(x, a) in O(1) using pi(x) lookup table if x < prime[a+1]^2.
///        * Cache results of phi(x, a) if x and a are small.
///        * Calculate all phi(x, a) = 1 upfront in O(1).
///        * Stop recursion at c instead of 1.
///
///       [1] Tomás Oliveira e Silva, Computing pi(x): the combinatorial
///           method, Revista do DETUA, vol. 4, no. 6, March 2006, p. 761.
///           http://sweet.ua.pt/tos/bib/5.4.pdf
///       [2] phi(x, a) = (x / pp) * φ(pp) + phi(x % pp, a)
///           with pp = 2 * 3 * ... * prime[a] 
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef GENERATE_PHI_HPP
#define GENERATE_PHI_HPP

#include <primecount-internal.hpp>
#include <fast_div.hpp>
#include <imath.hpp>
#include <min.hpp>
#include <PhiTiny.hpp>
#include <PiTable.hpp>
#include <pod_vector.hpp>

#include <stdint.h>
#include <array>
#include <limits>

namespace {

using namespace std;
using namespace primecount;

template <typename Primes>
class PhiCache
{
public:
  PhiCache(int64_t limit,
           const Primes& primes,
           const PiTable& pi)
    : primes_(primes),
      pi_(pi)
  {
    // We cache phi(x, a) results if x <= cache_limit_ (and a <= 100).
    // Actually we cache phi(x, a) results if (x + 1) / 2 <= cache_limit_
    // because phi(x, a) only changes its result if x is odd (for the
    // same a). This trick allows us to double the capacity of our cache
    // without increasing its memory usage.
    cache_limit_ = numeric_limits<cache_t>::max();
    cache_limit_ = min(cache_limit_, isqrt(limit));
  }

  /// Calculate phi(x, a) using the recursive formula:
  /// phi(x, a) = phi(x, a - 1) - phi(x / primes[a], a - 1)
  ///
  template <int SIGN>
  int64_t phi(int64_t x, int64_t a)
  {
    if (x <= (int64_t) primes_[a])
      return SIGN;
    else if (is_phi_tiny(a))
      return phi_tiny(x, a) * SIGN;
    else if (is_pix(x, a))
      return (pi_[x] - a + 1) * SIGN;
    else if (is_cached(x, a))
      return phi_cache(x, a) * SIGN;

    int64_t sqrtx = isqrt(x);
    int64_t pi_sqrtx = a;
    int64_t c = PhiTiny::get_c(sqrtx);
    int64_t sum = 0;

    if (sqrtx < pi_.size())
      pi_sqrtx = min(pi_[sqrtx], a);

    // Move out of the loop the calculations where phi(xp, i) = 1
    // phi(x, a) = 1 if primes[a] >= x
    // xp = x / primes[i + 1]
    // phi(xp, i) = 1 if primes[i] >= x / primes[i + 1]
    // phi(xp, i) = 1 if primes[i] >= sqrt(x)
    // phi(xp, i) = 1 if i >= pi(sqrt(x))
    // \sum_{i = pi(sqrt(x))}^{a - 1} phi(xp, i) = a - pi(sqrt(x))
    //
    sum += (pi_sqrtx - a) * SIGN;
    sum += phi_tiny(x, c) * SIGN;

    int64_t i;
    int64_t xp;

    for (i = c; i < pi_sqrtx; i++)
    {
      xp = fast_div(x, primes_[i + 1]);
      if (is_pix(xp, i))
        goto use_faster_pix;

      sum += phi<-SIGN>(xp, i);
    }

    for (; i < pi_sqrtx; i++)
    {
      xp = fast_div(x, primes_[i + 1]);
      use_faster_pix:;
      sum += (pi_[xp] - i + 1) * -SIGN;
    }

    update_cache(x, a, sum);

    return sum;
  }

private:
  /// Cache phi(x, a) results if (x + 1) / 2 <= cache_limit_
  uint64_t cache_limit_ = 0;
  enum { MAX_A = 100 };
  using cache_t = uint16_t;
  array<vector<cache_t>, MAX_A> cache_;
  const Primes& primes_;
  const PiTable& pi_;

  bool is_pix(int64_t x, int64_t a) const
  {
    return x < pi_.size() &&
           x < isquare(primes_[a + 1]);
  }

  bool is_cached(uint64_t x, uint64_t a) const
  {
    x = ceil_div(x, 2);
    return a < cache_.size() &&
           x < cache_[a].size() &&
           cache_[a][x] != 0;
  }

  int64_t phi_cache(uint64_t x, uint64_t a) const
  {
    // We cache phi(x, a) results if x <= cache_limit_ (and a <= 100).
    // Actually we cache phi(x, a) results if (x + 1) / 2 <= cache_limit_
    // because phi(x, a) only changes its result if x is odd (for the
    // same a). This trick allows us to double the capacity of our cache
    // without increasing its memory usage.
    x = ceil_div(x, 2);
    return cache_[a][x];
  }

  void update_cache(uint64_t x, uint64_t a, int64_t sum)
  {
    x = ceil_div(x, 2);

    if (a >= cache_.size() ||
        x > cache_limit_)
      return;

    if (x >= cache_[a].size())
    {
      cache_[a].reserve(cache_limit_ + 1);
      cache_[a].resize(x + 1, 0);
    }

    sum = abs(sum);
    assert(sum <= numeric_limits<cache_t>::max());
    cache_[a][x] = (cache_t) sum;
  }
};

/// Returns a vector with phi(x, i - 1) values such that
/// phi[i] = phi(x, i - 1) for 1 <= i <= a.
/// phi(x, a) counts the numbers <= x that are not
/// divisible by any of the first a primes.
///
template <typename Primes>
pod_vector<int64_t>
generate_phi(int64_t x,
             int64_t a,
             const Primes& primes,
             const PiTable& pi)
{
  int64_t size = a + 1;
  pod_vector<int64_t> phi(size);
  phi[0] = 0;

  if (size > 1)
  {
    if ((int64_t) primes[a] > x)
      a = pi[x];

    phi[1] = x;
    int64_t sqrtx = isqrt(x);
    int64_t pi_sqrtx = a;
    PhiCache<Primes> cache(x, primes, pi);

    if (sqrtx < pi.size())
      pi_sqrtx = min(pi[sqrtx] + 1, a);

    for (int64_t i = 2; i <= pi_sqrtx; i++)
      phi[i] = phi[i - 1] + cache.template phi<-1>(x / primes[i - 1], i - 2);

    for (int64_t i = pi_sqrtx + 1; i <= a; i++)
      phi[i] = phi[i - 1] - (x > 0);

    for (int64_t i = a + 1; i < size; i++)
      phi[i] = x > 0;
  }

  return phi;
}

} // namespace

#endif
