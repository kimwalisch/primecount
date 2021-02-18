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
#include <BitSieve240.hpp>
#include <generate.hpp>
#include <gourdon.hpp>
#include <fast_div.hpp>
#include <imath.hpp>
#include <min.hpp>
#include <PhiTiny.hpp>
#include <print.hpp>
#include <popcnt.hpp>

#include <stdint.h>
#include <cassert>
#include <array>
#include <utility>
#include <vector>

using namespace std;
using namespace primecount;

namespace {

class PhiCache : public BitSieve240
{
public:
  PhiCache(int64_t limit,
           const vector<int32_t>& primes,
           const PiTable& pi) :
    primes_(primes),
    pi_(pi)
  {
    // sieve_[a][x] uses at most max_x_bytes
    uint64_t max_x_bytes = 96 << 10;
    uint64_t numbers_per_sieve_byte = 240 / sizeof(uint64_t);
    max_x_ = max_x_bytes * numbers_per_sieve_byte;

    // This cache limit has been tuned for both pi_legendre(x) and
    // pi_meissel(x) as these functions are frequently used in primecount.
    // The idea is to allocate only a small amount of cache memory for
    // tiny computations (because initializing a lot of cache memory is
    // slow). On the other hand, for large and long running computations
    // we use the maximum amount of cache memory.
    uint64_t cache_limit = isqrt(limit);
    max_x_ = min(cache_limit, max_x_);
    max_x_size_ = ceil_div(max_x_, 240);

    // Make sure that there are no uninitialized
    // bits in the last sieve array item.
    max_x_ = max_x_size_ * 240 - 1;
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
           x / 240 < sieve_[a].size();
  }

  int64_t phi_cache(uint64_t x, uint64_t a) const
  {
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
      // correspond to the offsets { 1, 7, 11, 13, 17, 19, 23, 29}.
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
  enum { MAX_A = 100 };
  const uint64_t max_a_ = MAX_A;
  /// sieve_[a] contains only numbers (1 bits) that are
  /// not divisible by any of the first a primes.
  array<vector<uint64_t>, MAX_A + 1> sieve_;
  /// sieve_counts_[a][i] contains the count of numbers < i * 240 that
  /// are not divisible by any of the first a primes.
  array<vector<uint32_t>, MAX_A + 1> sieve_counts_;
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

/// pi(x) <= pix_upper(x)
/// pi(x) <= x / (log(x) - 1.1) + 5, for x >= 4.
/// We use x >= 10 and +10 as a safety buffer.
/// https://en.wikipedia.org/wiki/Prime-counting_function#Inequalities
///
int64_t pix_upper(int64_t x)
{
  if (x <= 10)
    return 4;

  double pix = x / (log((double) x) - 1.1);
  return (int64_t) pix + 10;
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
  if (a < 1) return x;

  // phi(x, a) = 1 if prime[a] >= x
  if (x > 0 && a > x / 2)
    return 1;

  if (is_phi_tiny(a))
    return phi_tiny(x, a);

  // phi(x, a) = 1 if a >= pi(x)
  if (a >= pix_upper(x))
    return 1;

  int64_t sqrtx = isqrt(x);

  // Fast (a > pi(sqrt(x)) check with decent accuracy
  if (a > pix_upper(sqrtx))
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
