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
#include <BitSieve128.hpp>
#include <fast_div.hpp>
#include <imath.hpp>
#include <min.hpp>
#include <PhiTiny.hpp>
#include <PiTable.hpp>
#include <pod_vector.hpp>

#include <stdint.h>
#include <array>
#include <vector>
#include <limits>

namespace {

using namespace std;
using namespace primecount;

template <typename Primes>
class PhiCache : public BitSieve128
{
public:
  PhiCache(int64_t limit,
           const Primes& primes,
           const PiTable& pi) :
    primes_(primes),
    pi_(pi)
  {
    // sieve_[a][x] uses at most max_x_bytes
    uint64_t max_x_bytes = 128 << 10;
    uint64_t numbers_per_sieve_byte = 16;
    max_x_ = max_x_bytes * numbers_per_sieve_byte;

    // This cache limit has been tuned for both pi_legendre(x) and
    // pi_meissel(x) as these functions are frequently used in primecount.
    // The idea is to allocate only a small amount of cache memory for
    // tiny computations (because initializing a lot of cache memory is
    // slow). On the other hand, for large and long running computations
    // we use the maximum amount of cache memory.
    uint64_t cache_limit = (uint64_t) pow((double) limit, 1 / 2.5);
    max_x_ = min(cache_limit, max_x_);
    max_x_size_ = ceil_div(max_x_, 128);

    // Make sure that there are no uninitialized
    // bits in the last sieve array item.
    max_x_ = max_x_size_ * 128 - 1;
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
    int64_t i = c;

    for (; i < a; i++)
    {
      // phi(x / prime[i+1], prime[i]) = 1 if prime[i] * prime[i+1] >= x.
      // However we can do slightly better:
      // If prime[i + 1] > sqrt(x) then x / prime[i + 1] < sqrt(x).
      // Hence x / prime[i + 1] <= sqrt(x) - 1.
      // phi(sqrt(x) - 1, i) = 1 as prime[i] is the largest
      // prime <= (sqrt(x) - 1) as prime[i + 1] > sqrt(x) and hence there
      // is no other prime number inside [prime[i] + 1, sqrt(x) - 1].
      if (primes_[i + 1] > sqrtx)
        break;
      int64_t xp = fast_div(x, primes_[i + 1]);
      if (is_pix(xp, i))
        break;
      sum += phi<-SIGN>(xp, i);
    }

    for (; i < a; i++)
    {
      if (primes_[i + 1] > sqrtx)
        break;
      int64_t xp = fast_div(x, primes_[i + 1]);
      sum += (pi_[xp] - i + 1) * -SIGN;
    }

    // phi(x, a) = 1 for all primes[a] >= x
    sum += (a - i) * -SIGN;
    sieve_cache(a);

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
    return a < sieve_.size() &&
           x / 128 < sieve_[a].size();
  }

  int64_t phi_cache(uint64_t x, uint64_t a) const
  {
    assert(x / 128 < sieve_[a].size());
    assert(x / 128 < sieve_counts_[a].size());
    uint64_t bitmask = unset_bits_[x % 128];
    uint64_t bits = sieve_[a][x / 128];
    return sieve_counts_[a][x / 128] + popcnt64(bits & bitmask);
  }

  /// Eratosthenes-like sieving algorithm that removes the first a primes
  /// and their multiples from the sieve array. Additionally this
  /// algorithm counts the numbers that are not divible by any of the
  /// first a primes after sieving has completed. The sieve array and the
  /// sieve_counts array later serve us as a phi(x, a) cache.
  ///
  void sieve_cache(uint64_t a)
  {
    if (a <= prev_a_ ||
        a >= sieve_.size())
      return;

    uint64_t i = prev_a_ + 1;
    uint64_t tiny_a = PhiTiny::max_a();
    prev_a_ = a;

    for (; i <= a; i++)
    {
      if (primes_[i] == 2)
        sieve_[i].resize(max_x_size_, ~0ull);
      else
      {
        // Initalize phi(x, a) with phi(x, a - 1) then
        // remove prime[a] and its multiples.
        sieve_[i] = sieve_[i - 1];
        uint64_t prime = primes_[i];

        if (prime <= max_x_)
          sieve_[i][prime / 128] &= unset_bit_[prime % 128];
        for (uint64_t n = prime * prime; n <= max_x_; n += prime * 2)
          sieve_[i][n / 128] &= unset_bit_[n % 128];
      }

      if (i > tiny_a)
      {
        uint64_t sum = 0;
        sieve_counts_[i].reserve(max_x_size_);
        for (uint64_t j = 0; j < max_x_size_; j++)
        {
          // Count of the numbers < j * 128 that are not
          // divisible by any of the first i primes.
          sieve_counts_[i].push_back(sum);
          sum += popcnt64(sieve_[i][j]);
        }
      }
    }
  }

  uint64_t max_x_ = 0;
  uint64_t max_x_size_ = 0;
  uint64_t prev_a_ = 0;
  enum { MAX_A = 100 };
  /// sieve_[a] contains only numbers (1 bits) that are
  /// not divisible by any of the first a primes.
  array<vector<uint64_t>, MAX_A> sieve_;
  /// sieve_counts_[a][i] contains the count of numbers < i * 128 that
  /// are not divisible by any of the first a primes.
  array<vector<uint32_t>, MAX_A> sieve_counts_;
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
    int64_t i = 2;
    int64_t sqrtx = isqrt(x);
    PhiCache<Primes> cache(x, primes, pi);

    // 2 <= i <= pi(sqrt(x)) + 1
    for (; i <= a && primes[i - 1] <= sqrtx; i++)
      phi[i] = phi[i - 1] + cache.template phi<-1>(x / primes[i - 1], i - 2);

    // pi(sqrt(x)) + 1 < i <= a
    for (; i <= a; i++)
      phi[i] = phi[i - 1] - (x > 0);

    // a < i < size
    for (; i < size; i++)
      phi[i] = x > 0;
  }

  return phi;
}

} // namespace

#endif
