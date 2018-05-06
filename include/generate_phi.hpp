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

#ifndef GENERATE_PHI_HPP
#define GENERATE_PHI_HPP

#include <primecount-internal.hpp>
#include <fast_div.hpp>
#include <imath.hpp>
#include <PhiTiny.hpp>
#include <PiTable.hpp>

#include <stdint.h>
#include <array>
#include <vector>
#include <limits>

namespace {

using namespace std;
using namespace primecount;

/// Cache phi(x, a) results if a < MAX_A
const int MAX_A = 100;

template <typename Primes>
class PhiCache
{
public:
  PhiCache(Primes& primes, PiTable& pi)
    : primes_(primes),
      pi_(pi)
  { }

  /// Calculate phi(x, a) using the recursive formula:
  /// phi(x, a) = phi(x, a - 1) - phi(x / prime(a), a - 1)
  ///
  template <int SIGN>
  int64_t phi(int64_t x, int64_t a)
  {
    if (x <= prime(a))
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

    // Move out of the loop the calculations where phi(x2, i) = 1
    // phi(x, a) = 1 if prime(a) >= x
    // x2 = x / prime(i + 1)
    // phi(x2, i) = 1 if prime(i) >= x / prime(i + 1)
    // phi(x2, i) = 1 if prime(i) >= sqrt(x)
    // phi(x2, i) = 1 if i >= pi(sqrt(x))
    // \sum_{i = pi(sqrt(x))}^{a - 1} phi(x2, i) = a - pi(sqrt(x))
    //
    sum += (pi_sqrtx - a) * SIGN;
    sum += phi_tiny(x, c) * SIGN;

    for (int64_t i = c; i < pi_sqrtx; i++)
    {
      int64_t x2 = fast_div(x, prime(i + 1));

      if (is_pix(x2, i))
        sum += (pi_[x2] - i + 1) * -SIGN;
      else
        sum += phi<-SIGN>(x2, i);
    }

    update_cache(x, a, sum);

    return sum;
  }

private:
  using T = uint16_t;
  array<vector<T>, MAX_A> cache_;
  Primes& primes_;
  PiTable& pi_;

  int64_t prime(int64_t i) const
  {
    return primes_[i];
  }

  void update_cache(uint64_t x, uint64_t a, int64_t sum)
  {
    if (a < cache_.size() &&
        x <= numeric_limits<T>::max())
    {
      if (x >= cache_[a].size())
        cache_[a].resize(x + 1, 0);

      cache_[a][x] = (T) abs(sum);
    }
  }

  bool is_pix(int64_t x, int64_t a) const
  {
    return x < pi_.size() &&
           x < isquare(prime(a + 1));
  }

  bool is_cached(uint64_t x, uint64_t a) const
  {
    return a < cache_.size() && 
           x < cache_[a].size() && 
           cache_[a][x];
  }
};

/// Returns a vector with phi(x, i - 1) values such that
/// phi[i] = phi(x, i - 1) for 1 <= i <= a.
/// phi(x, a) counts the numbers <= x that are not
/// divisible by any of the first a primes.
///
template <typename Primes>
vector<int64_t> generate_phi(int64_t x, int64_t a, Primes& primes, PiTable& pi)
{
  int64_t size = a + 1;

  if (primes[a] > x)
    a = pi[x];

  vector<int64_t> phi;
  phi.reserve(size);
  phi.resize(a + 1, (x > 0) * -1);
  phi.resize(size, x > 0);

  if (size > 1)
  {
    phi[1] = x;
    int64_t sqrtx = isqrt(x);
    int64_t pi_sqrtx = a;
    PhiCache<Primes> cache(primes, pi);

    if (sqrtx < pi.size())
      pi_sqrtx = min(pi[sqrtx] + 1, a);

    for (int64_t i = 2; i <= pi_sqrtx; i++)
      phi[i] = cache.template phi<-1>(x / primes[i - 1], i - 2);

    // calculate phi(x, a) using partial results
    for (int64_t i = 2; i <= a; i++)
      phi[i] += phi[i - 1];
  }

  return phi;
}

} // namespace

#endif
