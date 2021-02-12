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
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
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
#include <ModuloWheel.hpp>
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
    // The cache limit has been tuned for both pi_legendre(x) and
    // pi_meissel(x) as these functions are frequently used in primecount.
    // The idea is to allocate only a small amount of cache memory for
    // tiny computations (because initializing a lot of cache memory is
    // slow). On the other hand, for large and long running computations
    // we use the maximum amount of cache memory.
    auto u16_max = numeric_limits<uint16_t>::max();
    cache_limit_ = (uint64_t) pow((double) limit, 1 / 2.5);
    cache_limit_ = min(cache_limit_, u16_max);
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
    int64_t c = PhiTiny::get_c(sqrtx);
    int64_t sum = phi_tiny(x, c) * SIGN;
    int64_t i, xp;

    for (i = c; i < a; i++)
    {
      xp = fast_div(x, primes_[i + 1]);
      if (is_pix(xp, i) ||
          xp <= primes_[i])
        break;
      sum += phi<-SIGN>(xp, i);
    }

    for (; i < a; i++)
    {
      xp = fast_div(x, primes_[i + 1]);
      if (xp <= primes_[i])
        break;
      sum += (pi_[xp] - i + 1) * -SIGN;
    }

    // phi(x, a) = 1 for all primes[a] >= x
    sum += (a - i) * -SIGN;
    update_cache(x, a, sum);

    return sum;
  }

private:
  /// phi(x, a) counts the numbers <= x that are not divisible by any of
  /// the first a primes. If x < prime[a+1]^2 then phi(x, a) counts the
  /// number of primes <= x, minus the first a primes, plus the number 1.
  /// Hence if x < prime[a+1]^2: phi(x, a) = pi(x) - a + 1.
  ///
  bool is_pix(int64_t x, int64_t a) const
  {
    return x < pi_.size() &&
           x < isquare(primes_[a + 1]);
  }

  bool is_cached(uint64_t x, uint64_t a) const
  {
    x = ModuloWheel::to_index(x);
    return a < cache_.size() &&
           x < cache_[a].size() &&
           cache_[a][x] != 0;
  }

  int64_t phi_cache(uint64_t x, uint64_t a) const
  {
    x = ModuloWheel::to_index(x);
    return cache_[a][x];
  }

  void update_cache(uint64_t x, uint64_t a, int64_t sum)
  {
    // In order to increase the capacity of our phi cache we use a
    // modulo 2310 wheel that does not store any numbers that are
    // divisible by 2, 3, 5, 7 and 11. In a modulo 2310 wheel 480 array
    // items span an interval of size 2310. Hence by using a modulo
    // 2310 wheel we can increase our cache's capacity by 2310 / 480
    // = 4.8125 without increasing its memory usage.
    x = ModuloWheel::to_index(x);

    if (a >= cache_.size() ||
        x > cache_limit_)
      return;

    if (x >= cache_[a].size())
    {
      cache_[a].reserve(cache_limit_ + 1);
      cache_[a].resize(x + 1, 0);
    }

    sum = abs(sum);
    assert(sum <= numeric_limits<uint16_t>::max());
    assert(cache_[a][x] == 0);
    cache_[a][x] = (uint16_t) sum;
  }

  uint64_t cache_limit_ = 0;
  enum { MAX_A = 100 };
  array<vector<uint16_t>, MAX_A> cache_;
  const Primes& primes_;
  const PiTable& pi_;
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
