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
///        * Cache results of phi(x, a)
///        * Calculate phi(x, a) using formula [2] if a <= 6
///        * Calculate phi(x, a) using pi(x) lookup table
///        * Calculate all phi(x, a) = 1 upfront
///        * Stop recursion at c instead of 1
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
  PhiCache(vector<int32_t>& primes,
           PiTable& pi) :
    primes_(primes),
    pi_(pi)
  { }

  /// Calculate phi(x, a) using the recursive formula:
  /// phi(x, a) = phi(x, a - 1) - phi(x / primes_[a], a - 1)
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
      return cache_[a][x] * SIGN;

    int64_t sqrtx = isqrt(x);
    int64_t pi_sqrtx = a;
    int64_t c = PhiTiny::get_c(sqrtx);
    int64_t sum = 0;

    if (sqrtx < pi_.size())
      pi_sqrtx = min(pi_[sqrtx], a);

    // Move out of the loop the calculations where phi(xp, i) = 1
    // phi(x, a) = 1 if primes_[a] >= x
    // xp = x / primes_[i + 1]
    // phi(xp, i) = 1 if primes_[i] >= x / primes_[i + 1]
    // phi(xp, i) = 1 if primes_[i] >= sqrt(x)
    // phi(xp, i) = 1 if i >= pi(sqrt(x))
    // \sum_{i = pi(sqrt(x))}^{a - 1} phi(xp, i) = a - pi(sqrt(x))
    //
    sum += (pi_sqrtx - a) * SIGN;
    sum += phi_tiny(x, c) * SIGN;

    for (int64_t i = c; i < pi_sqrtx; i++)
    {
      int64_t xp = fast_div(x, primes_[i + 1]);

      if (is_pix(xp, i))
        sum += (pi_[xp] - i + 1) * -SIGN;
      else
        sum += phi<-SIGN>(xp, i);
    }

    update_cache(x, a, sum);

    return sum;
  }

private:
  /// Cache phi(x, a) results if a < MAX_A
  enum { MAX_A = 100 };
  using U16 = uint16_t;
  array<vector<U16>, MAX_A> cache_;
  vector<int32_t>& primes_;
  PiTable& pi_;

  void update_cache(uint64_t x, uint64_t a, int64_t sum)
  {
    auto max_x = numeric_limits<U16>::max();

    if (a < cache_.size() &&
        x <= max_x)
    {
      if (x >= cache_[a].size())
      {
        uint64_t max_size = max_x + 1;
        uint64_t size = min(x * 2, max_size);
        cache_[a].reserve(size);
        cache_[a].resize(x + 1, 0);
      }

      cache_[a][x] = (U16) abs(sum);
    }
  }

  bool is_pix(int64_t x, int64_t a) const
  {
    return x < pi_.size() &&
           x < isquare(primes_[a + 1]);
  }

  bool is_cached(uint64_t x, uint64_t a) const
  {
    return a < cache_.size() && 
           x < cache_[a].size() && 
           cache_[a][x];
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
      PiTable pi(max(sqrtx, primes[a]), threads);
      PhiCache cache(primes, pi);

      int64_t c = PhiTiny::get_c(sqrtx);
      int64_t pi_sqrtx = min(pi[sqrtx], a);
      int64_t thread_threshold = (int64_t) 1e10;
      threads = ideal_num_threads(threads, x, thread_threshold);

      sum = phi_tiny(x, c) - a + pi_sqrtx;

      #pragma omp parallel for num_threads(threads) schedule(dynamic, 16) firstprivate(cache) reduction(+: sum)
      for (int64_t i = c; i < pi_sqrtx; i++)
        sum += cache.phi<-1>(x / primes[i + 1], i);
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
