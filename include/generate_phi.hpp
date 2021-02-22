///
/// @file  generate_phi.hpp
/// @brief The PhiCache class calculates the partial sieve function
///        (a.k.a. Legendre-sum) using the recursive formula:
///        phi(x, a) = phi(x, a - 1) - phi(x / primes[a], a - 1).
///        phi(x, a) counts the numbers <= x that are not divisible
///        by any of the first a primes. The algorithm used is an
///        optimized version of the recursive algorithm described in
///        Tomás Oliveira e Silva's paper [1]. I have added 5
///        optimizations to my implementation which speed up the
///        computation by several orders of magnitude:
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
#include <BitSieve240.hpp>
#include <fast_div.hpp>
#include <imath.hpp>
#include <min.hpp>
#include <PhiTiny.hpp>
#include <PiTable.hpp>
#include <pod_vector.hpp>
#include <popcnt.hpp>

#include <stdint.h>
#include <cassert>
#include <utility>
#include <vector>

namespace {

using namespace std;
using namespace primecount;

template <typename Primes>
class PhiCache : public BitSieve240
{
public:
  PhiCache(uint64_t x,
           uint64_t a,
           const Primes& primes,
           const PiTable& pi) :
    primes_(primes),
    pi_(pi)
  {
    // Cache phi(n, a) if n <= sqrt(x) && a <= max_a.
    // The value max_a = 100 has been determined empirically
    // by running benchmarks. Using a smaller or larger
    // max_a with the same amount of memory (max_megabytes)
    // decreases the performance.
    uint64_t max_a = 100;
    uint64_t sqrtx = isqrt(x);
    uint64_t tiny_a = PhiTiny::max_a();

    // Make sure we cache only frequently used values
    a = a - min(a, 30);
    max_a = min(a, max_a);

    if (max_a <= tiny_a)
      return;

    // The cache (i.e. the sieve and sieve_counts arrays)
    // uses at most max_megabytes per thread.
    uint64_t max_megabytes = 16;
    uint64_t indexes = max_a - tiny_a;
    uint64_t max_bytes = max_megabytes << 20;
    uint64_t max_bytes_per_index = max_bytes / indexes;
    uint64_t numbers_per_sieve_byte = 240 / sizeof(uint64_t);
    max_x_ = ((max_bytes_per_index * 2) / 3) * numbers_per_sieve_byte;
    max_x_ = min(sqrtx, max_x_);
    max_x_size_ = ceil_div(max_x_, 240);

    if (max_x_size_ == 0)
      return;

    // Make sure that there are no uninitialized
    // bits in the last sieve array element.
    max_x_ = max_x_size_ * 240 - 1;
    max_a_ = max_a;
    sieve_.resize(max_a_ + 1);
    sieve_counts_.resize(max_a_ + 1);
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

    // Cache all small phi(x, i) results with:
    // x <= max_x && i <= min(a, max_a)
    sieve_cache(x, a);

    int64_t sqrtx = isqrt(x);
    int64_t c = PhiTiny::get_c(sqrtx);
    int64_t larger_c = min(a, max_a_cached_);
    int64_t sum, i;

    if (c >= larger_c ||
        !is_cached(x, larger_c))
      sum = phi_tiny(x, c) * SIGN;
    else
    {
      c = larger_c;
      assert(larger_c <= a);
      sum = phi_cache(x, c) * SIGN;
    }
  
    for (i = c; i < a; i++)
    {
      // phi(x / prime[i+1], i) = 1 if x / prime[i+1] <= prime[i].
      // However we can do slightly better:
      // If prime[i+1] > sqrt(x) and prime[i] <= sqrt(x) then
      // phi(x / prime[i+1], i) = 1 even if x / prime[i+1] > prime[i].
      // This works because in this case there is no other prime
      // inside the interval ]prime[i], x / prime[i+1]].
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
    return x <= max_x_ &&
           a <= max_a_cached_;
  }

  int64_t phi_cache(uint64_t x, uint64_t a) const
  {
    assert(a < sieve_.size());
    assert(x / 240 < sieve_[a].size());

    uint64_t bitmask = unset_larger_[x % 240];
    uint64_t bits = sieve_[a][x / 240];
    return sieve_counts_[a][x / 240] + popcnt64(bits & bitmask);
  }

  /// Cache phi(x, i) results with: x <= max_x && i <= min(a, max_a).
  /// Eratosthenes-like sieving algorithm that removes the first a primes
  /// and their multiples from the sieve array. Additionally this
  /// algorithm counts the numbers that are not divible by any of the
  /// first a primes after sieving has completed. The sieve array and the
  /// sieve_counts array later serve us as a phi(x, a) cache.
  ///
  void sieve_cache(uint64_t x, uint64_t a)
  {
    a = min(a, max_a_);

    if (x > max_x_ ||
        a <= max_a_cached_)
      return;

    uint64_t i = max_a_cached_ + 1;
    uint64_t tiny_a = PhiTiny::max_a();
    max_a_cached_ = a;
    i = max(i, 3);

    for (; i <= a; i++)
    {
      // Each bit in the sieve array corresponds to an integer that
      // is not divisible by 2, 3 and 5. The 8 bits of each byte
      // correspond to the offsets { 1, 7, 11, 13, 17, 19, 23, 29 }.
      if (i == 3)
        sieve_[i].resize(max_x_size_, ~0ull);
      else
      {
        // Initalize phi(x, i) with phi(x, i - 1)
        if (i - 1 <= tiny_a)
          sieve_[i] = std::move(sieve_[i - 1]);
        else
          sieve_[i] = sieve_[i - 1];

        // Remove prime[i] and its multiples
        uint64_t prime = primes_[i];
        if (prime <= max_x_)
          sieve_[i][prime / 240] &= unset_bit_[prime % 240];
        for (uint64_t n = prime * prime; n <= max_x_; n += prime * 2)
          sieve_[i][n / 240] &= unset_bit_[n % 240];

        if (i > tiny_a)
        {
          uint64_t count = 0;
          sieve_counts_[i].reserve(max_x_size_);

          // Fill an array with the cumulative 1 bit counts.
          // sieve_counts_[i][j] contains the count of numbers < j * 240
          // that are not divisible by any of the first i primes.
          for (uint64_t j = 0; j < max_x_size_; j++)
          {
            sieve_counts_[i].push_back((uint32_t) count);
            count += popcnt64(sieve_[i][j]);
          }
        }
      }
    }
  }

  uint64_t max_x_ = 0;
  uint64_t max_x_size_ = 0;
  uint64_t max_a_cached_ = 0;
  uint64_t max_a_ = 0;
  /// sieve_[a] contains only numbers (1 bits) that are
  /// not divisible by any of the first a primes.
  vector<vector<uint64_t>> sieve_;
  /// sieve_counts_[a][i] contains the count of numbers < i * 240 that
  /// are not divisible by any of the first a primes.
  vector<vector<uint32_t>> sieve_counts_;
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
    PhiCache<Primes> cache(x, a, primes, pi);

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
