///
/// @file  nth_prime_sieve.hpp
/// @brief In the nth prime algorithm we first count the number of
///        primes up to an nth prime approximation. Next, we generate
///        primes using a special segmented sieve of Eratosthenes
///        algorithm with low memmory usage to find the actual nth
///        prime (which is close to the nth prime approximation).
///
///        Since we need to generate prime numbers close to the nth
///        prime which could potentially be as large as 10^30, we
///        cannot use the traditoinal segmented sieve of Eratosthenes
///        due to its O(n^(1/2)) memory usage. Therefore our
///        implementation uses a segment size of O(n^(1/3)) which
///        slightly deteriorates the runtime complexity of our
///        segmented sieve of Eratosthenes implementation. However,
///        our nth prime approximation is off by less than n^(1/2)
///        and therefore the slightly worse runtime complexity
///        of our sieving algorithm does not deteriorate the overall
///        runtime complexity of our nth prime algorithm.
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "BitSieve240.hpp"
#include "print.hpp"

#include <primecount.hpp>
#include <primecount-config.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <ctz.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <macros.hpp>
#include <popcnt.hpp>
#include <Vector.hpp>

#include <algorithm>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace {

using namespace primecount;

template <typename T>
class NthPrimeSieve : public BitSieve240
{
public:
  T get_low() const
  {
    return low_;
  }

  uint64_t get_count() const
  {
    return count_;
  }

  /// Sieve interval [low, high]
  template <typename UT>
  void sieve(UT low, UT high)
  {
    if (high < 2)
      return;

    // NthPrimeSieve cannot generate the primes 2, 3 and 5.
    // In our bit sieve array the 8 bits of each byte correspond
    // to the offsets: { 1, 7, 11, 13, 17, 19, 23, 29 }.
    // Hence, our bit sieve array excludes 2, 3 and 5.
    if (low <= 5)
      throw primecount_error("NthPrimeSieve: low must > 5");

    UT old_low = low;
    if (low % 240)
      low -= low % 240;

    low_ = low;
    UT dist = (high - low) + 1;
    uint64_t size = (uint64_t) ceil_div(dist, 240);
    uint64_t sqrt_high = (uint64_t) isqrt(high);
    uint64_t prime;

    sieve_.resize(size);
    std::fill(sieve_.begin(), sieve_.end(), ~0ull);
    sieve_.front() &= unset_smaller_[old_low % 240];
    sieve_.back() &= unset_larger_[high % 240];
    primesieve::iterator iter(7, sqrt_high);

    if (sqrt_high < low)
    {
      while ((prime = iter.next_prime()) <= sqrt_high)
      {
        // Calculate first multiple > low
        // that is not divisible by 2.
        UT q = low / prime;
        UT n = prime * (q + 1 + (q & 1));
        ASSERT(n > prime);
        ASSERT(n % 2 != 0);
        uint64_t i = (uint64_t) (n - low);
        uint64_t i_max = (uint64_t) (high - low);

        // Cross-off multiples
        for (; i <= i_max; i += prime * 2)
          sieve_[i / 240] &= unset_bit_[i % 240];
      }
    }
    else
    {
      while ((prime = iter.next_prime()) <= sqrt_high)
      {
        // Calculate first multiple > low
        // that is not divisible by 2.
        UT q = low / prime;
        UT n = prime * (q + 1 + (q & 1));
        n = max(n, UT(prime) * prime);
        ASSERT(n > prime);
        ASSERT(n % 2 != 0);
        uint64_t i = (uint64_t) (n - low);
        uint64_t i_max = (uint64_t) (high - low);

        // Cross-off multiples
        for (; i <= i_max; i += prime * 2)
          sieve_[i / 240] &= unset_bit_[i % 240];
      }
    }

    uint64_t count = 0;

    // Count primes (1 bits)
    for (uint64_t bits : sieve_)
      count += popcnt64(bits);

    count_ = count;
  }

  T find_nth_prime(uint64_t n) const
  {
    ASSERT(n > 0);
    ASSERT(n <= count_);

    uint64_t count = 0;

    for (std::size_t i = 0; i < sieve_.size(); i++)
    {
      uint64_t bits = sieve_[i];
      uint64_t count_bits = popcnt64(bits);

      if_likely(count + count_bits < n)
        count += count_bits;
      else
      {
        for (; bits; bits &= bits - 1)
        {
          if (++count == n)
          {
            uint64_t bit_index = ctz64(bits);
            uint64_t bit_value = bit_values_[bit_index];
            T prime = low_ + i * 240 + bit_value;
            return prime;
          }
        }
      }
    }

    return 0;
  }

private:
  T low_ = 0;
  uint64_t count_ = 0;
  Vector<uint64_t> sieve_;
};

/// The aligned_vector class aligns each of its
/// elements on a new cache line in order to avoid
/// false sharing (cache trashing) when multiple
/// threads write to adjacent elements.
///
template <typename T>
class aligned_vector
{
  static_assert(sizeof(T) < MAX_CACHE_LINE_SIZE,
                "sizeof(T) must be < MAX_CACHE_LINE_SIZE");

public:
  aligned_vector(std::size_t size) : vect_(size) { }
  std::size_t size() const { return vect_.size(); }
  T& operator[](std::size_t pos) { return vect_[pos].val; }

private:
  struct CacheLine {
    T val;
    // We cannot use alignas(MAX_CACHE_LINE_SIZE) for
    // the CacheLine struct as GCC does not yet support
    // alignas(n) with n > 128. Also alignas(n) for
    // over-aligned data and dynamic memory allocation
    // is only supported since C++17.
    MAYBE_UNUSED char pad[MAX_CACHE_LINE_SIZE - sizeof(T)];
  };

  Vector<CacheLine> vect_;
};

/// Find the nth prime using a prime sieve.
/// @sieve_forward = true:  Find nth prime >= nth_prime_approx.
/// @sieve_forward = false: Find nth prime <= nth_prime_approx.
///
template <bool sieve_forward, typename T>
T nth_prime_sieve(uint64_t n,
                  T nth_prime_approx,
                  int threads)
{
  ASSERT(n > 0);

  T nth_prime = 0;
  uint64_t count = 0;
  uint64_t while_iters = 0;

  T root3 = iroot<3>(nth_prime_approx);
  uint64_t thread_dist = (uint64_t) (root3 * 30);
  uint64_t min_thread_dist = 8 * 240;
  thread_dist = max(min_thread_dist, thread_dist);
  uint64_t avg_prime_gap = ilog(nth_prime_approx) + 2;
  uint64_t dist_approx = n * avg_prime_gap;

  threads = ideal_num_threads(dist_approx, threads, thread_dist);
  aligned_vector<NthPrimeSieve<T>> sieves(threads);
  bool print_vars = is_print();
  bool finished = false;
  double time;

  if (print_vars)
  {
    print("");
    print("=== nth_prime_sieve ===");
    print_nth_prime_sieve(n, sieve_forward, nth_prime_approx, dist_approx, thread_dist, threads);
    time = get_time();
  }

  #pragma omp parallel num_threads(threads)
  while (!finished)
  {
  #ifdef _OPENMP
    int thread_id = omp_get_thread_num();
  #else
    int thread_id = 0;
  #endif

    // Unsigned integer division is usually
    // faster than signed integer division.
    using UT = typename pstd::make_unsigned<T>::type;
    uint64_t i = while_iters * threads + thread_id;
    UT low = 0, high = 0;

    if (sieve_forward)
    {
      low = nth_prime_approx + i * thread_dist;
      high = low + thread_dist - 1;
    }
    else if ((UT) nth_prime_approx > i * thread_dist)
    {
      high = nth_prime_approx - i * thread_dist;
      low = (high - min(high, thread_dist)) + 1;
    }

    // Sieve the current segment [low, high].
    // If possible use fast 64-bit bit integer division
    // instead of slow 128-bit integer division.
    if ( low <= pstd::numeric_limits<uint64_t>::max() &&
        high <= pstd::numeric_limits<uint64_t>::max())
      sieves[thread_id].sieve((uint64_t) low, (uint64_t) high);
    else
      sieves[thread_id].sieve(low, high);

    // Wait until all threads have finished
    // computing their current segment.
    #pragma omp barrier
    #pragma omp single
    {
      while_iters++;

      for (int j = 0; j < threads; j++)
      {
        if (sieve_forward)
        {
          if (count + sieves[j].get_count() < n)
            count += sieves[j].get_count();
          else
          {
            // Nth prime is in the current segment
            nth_prime = sieves[j].find_nth_prime(n - count);
            finished = true;
            break;
          }
        }
        else // Sieve backwards
        {
          count += sieves[j].get_count();

          if (count >= n)
          {
            // Nth prime is in the current segment
            nth_prime = sieves[j].find_nth_prime((count - n) + 1);
            finished = true;
            break;
          }
          else if (sieves[j].get_low() == 0)
          {
            finished = true;
            break;
          }
        }
      }
    }
  }

  if (!nth_prime)
    throw primecount_error("Failed to find nth prime!");

  if (print_vars)
  {
    print("Status: 100%");
    print_seconds(get_time() - time);
  }

  return nth_prime;
}

} // namespace
