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
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef NTH_PRIME_SIEVE_HPP
#define NTH_PRIME_SIEVE_HPP

#include "BitSieve240.hpp"
#include "print.hpp"

#include <primecount.hpp>
#include <primecount-config.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <ctz.hpp>
#include <fast_div.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <macros.hpp>
#include <popcnt.hpp>
#include <Vector.hpp>

#include <algorithm>
#include <atomic>
#include <memory>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace {

using namespace primecount;

constexpr uint64_t min_iter_dist = uint64_t(1e8);

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
  void sieve(UT low, UT high, int threads)
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
    uint64_t dist = uint64_t((high - low) + 1);
    uint64_t size = ceil_div(dist, 240);

    if (size > sieve_size_)
    {
      sieve_size_ = size;
      sieve_.reset(new std::atomic<uint64_t>[sieve_size_]);
    }

    // Initialize sieve array, set all bits to 1
    auto* sieve = sieve_.get();
    for (uint64_t i = 0; i < sieve_size_; i++)
      sieve[i].store(~0ull, std::memory_order_relaxed);
    sieve[0].fetch_and(unset_smaller_[old_low % 240], std::memory_order_relaxed);
    sieve[sieve_size_ - 1].fetch_and(unset_larger_[high % 240], std::memory_order_relaxed);

    #pragma omp taskgroup
    {
      uint64_t sqrt_high = (uint64_t) isqrt(high);
      threads = ideal_num_threads(sqrt_high, threads, min_iter_dist);
      uint64_t thread_dist = ceil_div(sqrt_high, threads);

      for (int t = 1; t < threads; t++)
      {
        UT iter_start = thread_dist * t + 1;
        UT iter_stop = thread_dist * (t + 1);
        iter_start = max(iter_start, 7);
        iter_stop = min(iter_stop, sqrt_high);

        #pragma omp task
        cross_off(low, high, iter_start, iter_stop);
      }

      cross_off(low, high, 7, thread_dist);
    }

    uint64_t count = 0;

    // Count primes (1 bits)
    for (uint64_t i = 0; i < sieve_size_; i++)
      count += popcnt64(sieve[i].load(std::memory_order_relaxed));

    count_ = count;
  }

  /// Sieve interval [low, high]
  template <typename UT>
  void cross_off(UT low,
                 UT high,
                 uint64_t iter_start,
                 uint64_t iter_stop)
  {
    primesieve::iterator iter(iter_start, iter_stop);
    uint64_t i_max = uint64_t(high - low);
    uint64_t prime;
    auto* sieve = sieve_.get();

    if (isqrt(high) < low)
    {
      while ((prime = iter.next_prime()) <= iter_stop)
      {
        // Calculate first multiple > low
        // that is not divisible by 2.
        UT q = fast_div(low, prime);
        UT n = prime * (q + 1 + (q & 1));
        ASSERT(n > prime);
        ASSERT(n % 2 != 0);
        uint64_t i = uint64_t(n - low);

        // Cross-off multiples
        for (; i <= i_max; i += prime * 2)
          sieve[i / 240].fetch_and(unset_bit_[i % 240], std::memory_order_relaxed);
      }
    }
    else
    {
      while ((prime = iter.next_prime()) <= iter_stop)
      {
        // Calculate first multiple > low
        // that is not divisible by 2.
        UT q = fast_div(low, prime);
        UT n = prime * (q + 1 + (q & 1));
        n = max(n, UT(prime) * prime);
        ASSERT(n > prime);
        ASSERT(n % 2 != 0);
        uint64_t i = uint64_t(n - low);

        // Cross-off multiples
        for (; i <= i_max; i += prime * 2)
          sieve[i / 240].fetch_and(unset_bit_[i % 240], std::memory_order_relaxed);
      }
    }
  }

  T find_nth_prime(uint64_t n) const
  {
    ASSERT(n > 0);
    ASSERT(n <= count_);

    uint64_t count = 0;
    auto* sieve = sieve_.get();

    for (uint64_t i = 0; i < sieve_size_; i++)
    {
      uint64_t bits = sieve[i].load(std::memory_order_relaxed);
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
  uint64_t sieve_size_ = 0;
  std::unique_ptr<std::atomic<uint64_t>[]> sieve_;
};

/// The aligned_vector class aligns each of its
/// elements on a new cache line in order to avoid
/// false sharing (cache trashing) when multiple
/// threads write to adjacent elements.
///
template <typename T>
class aligned_vector
{
public:
  aligned_vector(std::size_t size) : vect_(size) { }
  std::size_t size() const { return vect_.size(); }
  T& operator[](std::size_t pos) { return vect_[pos].val; }

private:
  struct CacheLine {
    T val;
    // Round up the size of each element to a cache line
    // multiple so adjacent elements do not share cache lines.
    MAYBE_UNUSED char pad[MAX_CACHE_LINE_SIZE -
                          (sizeof(T) % MAX_CACHE_LINE_SIZE)];
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
                  int max_threads)
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

  int main_threads = ideal_num_threads(dist_approx, max_threads, thread_dist);
  int threads_per_segment = 1;

#if defined(_OPENMP)
  uint64_t sqrt_n = (uint64_t) isqrt(nth_prime_approx);
  int max_threads_per_segment = in_between(1, max_threads / main_threads, 16);
  threads_per_segment = ideal_num_threads(sqrt_n, max_threads_per_segment, min_iter_dist);
#endif

  aligned_vector<NthPrimeSieve<T>> sieves(main_threads);
  bool print_vars = is_print();
  bool finished = false;
  double time;

  if (print_vars)
  {
    print("");
    print("=== nth_prime_sieve ===");
    print_nth_prime_sieve(n, sieve_forward, nth_prime_approx, dist_approx,
        thread_dist, main_threads, threads_per_segment);
    time = get_time();
  }

  int total_threads = main_threads * threads_per_segment;
  #pragma omp parallel num_threads(total_threads)
  #pragma omp single
  {
    while (!finished)
    {
      uint64_t current_iter = while_iters++;

      #pragma omp taskgroup
      {
        for (int t = 0; t < main_threads; t++)
        {
          // Unsigned integer division is usually
          // faster than signed integer division.
          using UT = typename pstd::make_unsigned<T>::type;
          uint64_t i = current_iter * main_threads + t;
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

          #pragma omp task
          {
            // Sieve the current segment [low, high].
            // If possible use fast 64-bit bit integer division
            // instead of slow 128-bit integer division.
            if ( low <= pstd::numeric_limits<uint64_t>::max() &&
                high <= pstd::numeric_limits<uint64_t>::max())
              sieves[t].sieve((uint64_t) low, (uint64_t) high, threads_per_segment);
            else
              sieves[t].sieve(low, high, threads_per_segment);
          }
        }
      }

      for (int t = 0; t < main_threads; t++)
      {
        if (sieve_forward)
        {
          if (count + sieves[t].get_count() < n)
            count += sieves[t].get_count();
          else
          {
            // Nth prime is in the current segment
            nth_prime = sieves[t].find_nth_prime(n - count);
            finished = true;
            break;
          }
        }
        else // Sieve backwards
        {
          count += sieves[t].get_count();

          if (count >= n)
          {
            // Nth prime is in the current segment
            nth_prime = sieves[t].find_nth_prime((count - n) + 1);
            finished = true;
            break;
          }
          else if (sieves[t].get_low() == 0)
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

#endif
