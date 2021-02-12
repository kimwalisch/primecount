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
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <primecount-internal.hpp>
#include <generate.hpp>
#include <gourdon.hpp>
#include <fast_div.hpp>
#include <imath.hpp>
#include <min.hpp>
#include <ModuloWheel.hpp>
#include <PhiTiny.hpp>
#include <print.hpp>

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
           const vector<int32_t>& primes,
           const PiTable& pi) :
    primes_(primes),
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
    if (x <= primes_[a])
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
  const vector<int32_t>& primes_;
  const PiTable& pi_;
};

/// If a is very large (i.e. prime[a] > sqrt(x)) then we need to
/// calculate phi(x, a) using an alternative algorithm. First, because
/// in this case there actually exists a much faster algorithm. And
/// secondly, because storing the first a primes in a vector may use a
/// huge amount of memory and cause an out of memory error.
///
/// This alternative algorithm works if a >= pi(sqrt(x)). However, we
/// need to be very careful: phi_pix(x, a) may call pi_legendre(x) which
/// calls phi(x, a) with a = pi(sqrt(x)), which would then again call
/// phi_pix(x, a) thereby causing infinite recursion. In order to prevent
/// this issue this function must only be called with a > pi(sqrt(x)).
///
int64_t phi_pix(int64_t x, int64_t a, int threads)
{
  bool print = is_print();
  set_print(false);

  int64_t pix = pi(x, threads);
  int64_t sum;

  if (a <= pix)
    sum = pix - a + 1;
  else
    sum = 1;

  set_print(print);
  return sum;
}

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

  if (is_phi_tiny(a))
    return phi_tiny(x, a);

  int64_t sqrtx = isqrt(x);

  // Inaccurate but fast (a > pi(sqrt(x)) check
  if (a > sqrtx)
    return phi_pix(x, a, threads);

  PiTable pi(sqrtx, threads);
  int64_t pi_sqrtx = pi[sqrtx];

  // We use if (a > pi(sqrt(x)) here instead of (a >= pi(sqrt(x)) because
  // we want to prevent that our pi_legendre(x) uses this code path.
  // Otherwise pi_legendre(x) would switch to using pi_gourdon(x) under
  // the hood which is not what users expect. Also using (a >= pi(sqrt(x))
  // here would cause infinite recursion, more info at phi_pix(x, a).
  if (a > pi_sqrtx)
    return phi_pix(x, a, threads);

  auto primes = generate_n_primes<int32_t>(a);
  int64_t c = PhiTiny::get_c(sqrtx);
  int64_t sum = phi_tiny(x, c);
  int64_t thread_threshold = (int64_t) 1e10;
  threads = ideal_num_threads(threads, x, thread_threshold);

  #pragma omp parallel num_threads(threads) reduction(+: sum)
  {
    // Each thread uses its own PhiCache object in
    // order to avoid thread synchronization.
    PhiCache cache(x, primes, pi);

    #pragma omp for nowait schedule(dynamic, 16)
    for (int64_t i = c; i < a; i++)
      sum += cache.phi<-1>(x / primes[i + 1], i);
  }

  return sum;
}

/// The default phi(x, a) implementation does not print anything to the
/// screen as it is used by pi_simple(x) which is used all over the
/// place (e.g. to initialize S1, S2, P2, P3, ...) and we don't want to
/// print any info about this. Hence we also provide phi_print(x, a) for
/// use cases where we do want to print the result of phi(x, a).
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
