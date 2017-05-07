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
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <primecount-internal.hpp>
#include <generate.hpp>
#include <imath.hpp>
#include <PhiTiny.hpp>
#include <fast_div.hpp>
#include <min_max.hpp>

#include <stdint.h>
#include <algorithm>
#include <array>
#include <vector>
#include <limits>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;

namespace {

/// Cache phi(x, a) results if a < MAX_A
const int MAX_A = 100;

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
    else if (is_cached(x, a))
      return cache_[a][x] * SIGN;
    else if (is_pix(x, a))
      return (pi_[x] - a + 1) * SIGN;

    int64_t sqrtx = isqrt(x);
    int64_t pi_sqrtx = a;
    int64_t c = PhiTiny::get_c(sqrtx);
    int64_t sum = 0;

    if (sqrtx < pi_.size() && sqrtx < primes_[a])
      pi_sqrtx = pi_[sqrtx];

    // Move out of the loop the calculations where phi(x2, i) = 1
    // phi(x, a) = 1 if primes_[a] >= x
    // x2 = x / primes_[i + 1]
    // phi(x2, i) = 1 if primes_[i] >= x / primes_[i + 1]
    // phi(x2, i) = 1 if primes_[i] >= sqrt(x)
    // phi(x2, i) = 1 if i >= pi(sqrt(x))
    // \sum_{i = pi(sqrt(x))}^{a - 1} phi(x2, i) = a - pi(sqrt(x))
    //
    sum += (pi_sqrtx - a) * SIGN;
    sum += phi_tiny(x, c) * SIGN;

    for (int64_t i = c; i < pi_sqrtx; i++)
    {
      int64_t x2 = fast_div(x, primes_[i + 1]);

      if (is_pix(x2, i))
        sum += (pi_[x2] - i + 1) * -SIGN;
      else
        sum += phi<-SIGN>(x2, i);
    }

    update_cache(x, a, sum * SIGN);

    return sum;
  }

private:
  using cache_t = vector<uint16_t>;
  array<cache_t, MAX_A> cache_;
  vector<int32_t>& primes_;
  PiTable& pi_;

  void update_cache(uint64_t x, uint64_t a, uint64_t sum)
  {
    if (a < cache_.size() &&
        x <= numeric_limits<uint16_t>::max())
    {
      if (x >= cache_[a].size())
        cache_[a].resize(x + 1, 0);

      cache_[a][x] = (uint16_t) sum;
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

  print("");
  print("=== phi(x, a) ===");
  print("Count the numbers <= x coprime to the first a primes");

  double time = get_wtime();
  int64_t sum = 0;

  if (is_phi_tiny(a))
    sum = phi_tiny(x, a);
  else
  {
    auto primes = generate_n_primes<int32_t>(a);

    if (primes.at(a) >= x)
      sum = 1;
    else
    {
      // use a large pi(x) lookup table for speed
      int64_t sqrtx = isqrt(x);
      PiTable pi(max(sqrtx, primes[a]));
      PhiCache cache(primes, pi);

      int64_t pi_sqrtx = min(pi[sqrtx], a); 
      sum = x - a + pi_sqrtx;

      int64_t thread_threshold = ipow(10ll, 10);
      threads = ideal_num_threads(threads, x, thread_threshold);

      #pragma omp parallel for schedule(dynamic, 16) \
          num_threads(threads) firstprivate(cache) reduction(+: sum)
      for (int64_t i = 0; i < pi_sqrtx; i++)
        sum += cache.phi<-1>(x / primes[i + 1], i);
    }
  }

  print("phi", sum, time);
  return sum;
}

} // namespace
