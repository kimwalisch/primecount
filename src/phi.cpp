///
/// @file  phi.cpp
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

#include <PiTable.hpp>
#include <primecount-internal.hpp>
#include <print.hpp>
#include <generate.hpp>
#include <imath.hpp>
#include <PhiTiny.hpp>
#include <fast_div.hpp>
#include <min.hpp>

#include <stdint.h>
#include <algorithm>
#include <array>
#include <vector>
#include <limits>

using namespace std;
using namespace primecount;

namespace {

class PhiCache
{
public:
  PhiCache(int64_t limit,
           vector<int32_t>& primes,
           PiTable& pi) :
    primes_(primes),
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
    if (x <= primes_[a])
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
  vector<int32_t>& primes_;
  PiTable& pi_;

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

} // namespace

namespace primecount {

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
///
int64_t phi(int64_t x, int64_t a, int threads)
{
  if (x < 1) return 0;
  if (a > x) return 1;
  if (a < 1) return x;

  int64_t sum = 0;

  if (is_phi_tiny(a))
    sum = phi_tiny(x, a);
  else
  {
    auto primes = generate_n_primes<int32_t>(a);

    if (primes[a] >= x)
      sum = 1;
    else
    {
      // use large pi(x) lookup table for speed
      int64_t sqrtx = isqrt(x);
      int64_t c = PhiTiny::get_c(sqrtx);
      PiTable pi(max(sqrtx, primes[a]), threads);
      int64_t pi_sqrtx = min(pi[sqrtx], a);
      int64_t thread_threshold = (int64_t) 1e10;
      threads = ideal_num_threads(threads, x, thread_threshold);

      sum = phi_tiny(x, c) - a + pi_sqrtx;

      #pragma omp parallel num_threads(threads) reduction(+: sum)
      {
        // Each thread uses its own PhiCache object in
        // order to avoid thread synchronization.
        PhiCache cache(x, primes, pi);

        #pragma omp for nowait schedule(dynamic, 16)
        for (int64_t i = c; i < pi_sqrtx; i++)
          sum += cache.phi<-1>(x / primes[i + 1], i);
      }
    }
  }

  return sum;
}

/// The default phi(x, a) implementation does not print anything
/// to the screen as it is used by pi_simple(x) which is used
/// all over the place (e.g. to initialize S1, S2, P2, P3, ...)
/// and we don't want to print any info about this.
/// Hence we also provide phi_print(x, a) for use cases where we
/// do want to print the result of phi(x, a).
///
int64_t phi_print(int64_t x, int64_t a, int threads)
{
  print("");
  print("=== phi(x, a) ===");

  double time = get_time();
  int64_t sum = phi(x, a, threads);
  print("phi", sum, time);

  return sum;
}

} // namespace
