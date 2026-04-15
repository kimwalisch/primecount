///
/// @file  nth_prime_sieve.cpp
/// @brief In the nth prime algorithm we first count the number of
///        primes up to an nth prime approximation. Next, we generate
///        primes using a special segmented sieve of Eratosthenes
///        algorithm with low memory usage to find the actual nth
///        prime (which is close to the nth prime approximation).
///
///        Since we need to generate prime numbers close to the nth
///        prime which could potentially be as large as 10^30, we
///        cannot use the traditional segmented sieve of Eratosthenes
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

#include "nth_prime_sieve.hpp"
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
#include <RelaxedAtomic.hpp>
#include <Vector.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <memory>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace {

using namespace primecount;

/// NthPrimeSieve1 uses multi-threading with 1 thread per segment
/// whereas NthPrimeSieve2 uses multiple threads per segment.
/// NthPrimeSieve1 runs faster than NthPrimeSieve2 if the number
/// of segments is greater than or equal to the number of threads,
/// hence if threads_per_segment = 1.
///
template <typename T>
class NthPrimeSieve1 : public BitSieve240
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

    // NthPrimeSieve1 cannot generate the primes 2, 3 and 5.
    // In our bit sieve array the 8 bits of each byte correspond
    // to the offsets: { 1, 7, 11, 13, 17, 19, 23, 29 }.
    // Hence, our bit sieve array excludes 2, 3 and 5.
    if (low <= 5)
      throw primecount_error("NthPrimeSieve1: low must > 5");

    UT old_low = low;
    if (low % 240)
      low -= low % 240;

    low_ = low;
    uint64_t dist = uint64_t((high - low) + 1);
    uint64_t size = ceil_div(dist, 240);
    sieve_.resize(size);

    // Initialize sieve array, set all bits to 1
    std::fill(sieve_.begin(), sieve_.end(), ~0ull);
    sieve_.front() &= unset_smaller_[old_low % 240];
    sieve_.back() &= unset_larger_[high % 240];

    uint64_t prime;
    uint64_t i_max = uint64_t(high - low);
    uint64_t sqrt_high = (uint64_t) isqrt(high);
    primesieve::iterator iter(7, sqrt_high);

    if (sqrt_high < low)
    {
      while ((prime = iter.next_prime()) <= sqrt_high)
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
          sieve_[i / 240] &= unset_bit_[i % 240];
      }
    }
    else
    {
      while ((prime = iter.next_prime()) <= sqrt_high)
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
/// false sharing (cache thrashing) when multiple
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
T nth_prime_sieve1(uint64_t n,
                   T nth_prime_approx,
                   int threads)
{
  ASSERT(n > 0);
  ASSERT(nth_prime_approx > 0);

  // Calculate dist_approx such that:
  // nth_prime < nth_prime_approx + dist_approx
  double logx = std::log((double) nth_prime_approx);
  double log10_n = std::log10(max(n, 10));
  uint64_t dist_approx = (uint64_t) std::ceil(n * (logx + 1));
  dist_approx += uint64_t(logx * logx * log10_n);

  uint64_t root3 = (uint64_t) iroot<3>(nth_prime_approx);
  uint64_t max_thread_dist = uint64_t(root3 * 30);
  uint64_t thread_dist = in_between(240u, dist_approx, max_thread_dist);
  threads = ideal_num_threads(dist_approx, threads, thread_dist);
  aligned_vector<NthPrimeSieve1<T>> sieves(threads);
  double time;

  if (is_print())
  {
    print("");
    print("=== nth_prime_sieve ===");
    print_nth_prime_sieve(n, sieve_forward, nth_prime_approx,
        dist_approx, thread_dist, threads);
    time = get_time();
  }

  T nth_prime = 0;
  uint64_t count = 0;
  uint64_t while_iters = 0;
  bool finished = false;

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
    // If possible use fast 64-bit integer division
    // instead of slow 128-bit integer division.
    if ( low <= pstd::numeric_limits<uint64_t>::max() &&
        high <= pstd::numeric_limits<uint64_t>::max())
      sieves[thread_id].sieve(uint64_t(low), uint64_t(high));
    else
      sieves[thread_id].sieve(low, high);

    // Wait until all threads have finished
    // computing their current segment.
    #pragma omp barrier
    #pragma omp single
    {
      while_iters++;

      for (int t = 0; t < threads; t++)
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

  if_unlikely(!nth_prime)
    throw primecount_error("Failed to find nth prime!");

  if (is_print())
    print_seconds(get_time() - time);

  return nth_prime;
}

/// nth_prime_sieve2() requires OpenMP 4.0 or later
#if _OPENMP >= 201307

/// This thread threshold is used to calculate the number of
/// threads used in NthPrimeSieve2. Using a smaller thread
/// threshold typically increases the number of threads. For
/// small computations < 10^12 for which only a single thread
/// would be used by NthPrimeSieve2, we switch to the faster
/// NthPrimeSieve1 without atomic memory accesses.
///
constexpr uint64_t thread_threshold = uint64_t(2e7);

/// NthPrimeSieve2 is virtually identical to NthPrimeSieve1
/// except that NthPrimeSieve2 uses multiple threads per segment
/// whereas NthPrimeSieve1 uses a single thread per segment.
/// Because NthPrimeSieve2 uses multiple threads per segment, it
/// uses atomic memory accesses when crossing off multiples in
/// the sieve array. NthPrimeSieve2 is faster than NthPrimeSieve1
/// if the total number of segments is smaller than the number
/// of threads used.
///
template <typename T>
class NthPrimeSieve2 : public BitSieve240
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

    // NthPrimeSieve2 cannot generate the primes 2, 3 and 5.
    // In our bit sieve array the 8 bits of each byte correspond
    // to the offsets: { 1, 7, 11, 13, 17, 19, 23, 29 }.
    // Hence, our bit sieve array excludes 2, 3 and 5.
    if (low <= 5)
      throw primecount_error("NthPrimeSieve2: low must > 5");

    UT old_low = low;
    if (low % 240)
      low -= low % 240;

    low_ = low;
    uint64_t dist = uint64_t((high - low) + 1);
    sieve_size_ = ceil_div(dist, 240);

    if (sieve_size_ > sieve_capacity_)
    {
      sieve_.reset(new std::atomic<uint64_t>[sieve_size_]);
      sieve_capacity_ = sieve_size_;
    }

    // Initialize sieve array, set all bits to 1
    auto* sieve = sieve_.get();
    for (uint64_t i = 0; i < sieve_size_; i++)
      sieve[i].store(~0ull, std::memory_order_relaxed);
    sieve[0].fetch_and(unset_smaller_[old_low % 240], std::memory_order_relaxed);
    sieve[sieve_size_ - 1].fetch_and(unset_larger_[high % 240], std::memory_order_relaxed);

    #pragma omp taskgroup
    {
      // Calculate the number of threads, the thread chunk
      // size and the total number of thread chunks. These
      // values are important for load balancing, the 
      // number of chunks should be larger than the number
      // of threads to prevent load imbalance.
      uint64_t chunk_count = 1;
      uint64_t sqrt_high = (uint64_t) isqrt(high);
      uint64_t high_14 = isqrt(sqrt_high);
      uint64_t min_chunk_dist = max(high_14, uint64_t(1e6));
      threads = ideal_num_threads(sqrt_high, threads, thread_threshold);

      if (threads > 1)
      {
        chunk_count = sqrt_high / min_chunk_dist;
        double log10_chunk_count = std::log10(max(chunk_count, 10));
        double thread_chunks = std::pow(log10_chunk_count, 3);
        uint64_t max_chunk_count = int(thread_chunks) * threads;
        chunk_count = in_between(1u, chunk_count, max_chunk_count);
      }

      RelaxedAtomic<uint64_t> next_chunk(0);
      uint64_t chunk_dist = sqrt_high / chunk_count;
      chunk_dist = max(min_chunk_dist, chunk_dist);
      chunk_count = ceil_div(sqrt_high, chunk_dist);
      threads = (int) min(chunk_count, threads);

      // The main thread starts the worker threads
      for (int t = 1; t < threads; t++)
      {
        #pragma omp task shared(next_chunk)
        for (uint64_t i = next_chunk++; i < chunk_count; i = next_chunk++)
        {
          uint64_t iter_start = chunk_dist * i + 1;
          uint64_t iter_stop = chunk_dist * (i + 1);
          iter_start = max(iter_start, 7);
          iter_stop = min(iter_stop, sqrt_high);

          if (iter_start <= iter_stop)
            cross_off(low, high, iter_start, iter_stop);
        }
      }

      // After having started the worker threads the main
      // thread acts as an additional worker thread.
      for (uint64_t i = next_chunk++; i < chunk_count; i = next_chunk++)
      {
        uint64_t iter_start = chunk_dist * i + 1;
        uint64_t iter_stop = chunk_dist * (i + 1);
        iter_start = max(iter_start, 7);
        iter_stop = min(iter_stop, sqrt_high);

        if (iter_start <= iter_stop)
          cross_off(low, high, iter_start, iter_stop);
      }
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
    uint64_t sqrt_high = (uint64_t) isqrt(high);
    uint64_t i_max = uint64_t(high - low);
    uint64_t prime;

    auto* sieve = sieve_.get();

    if (sqrt_high < low)
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
        {
          uint64_t bitmask = unset_bit_[i % 240];
          if (bitmask != uint64_t(~0ull))
            sieve[i / 240].fetch_and(bitmask, std::memory_order_relaxed);
        }
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
        {
          uint64_t bitmask = unset_bit_[i % 240];
          if (bitmask != uint64_t(~0ull))
            sieve[i / 240].fetch_and(bitmask, std::memory_order_relaxed);
        }
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
  uint64_t sieve_capacity_ = 0;
  std::unique_ptr<std::atomic<uint64_t>[]> sieve_;
};

/// Find the nth prime using a prime sieve.
/// @sieve_forward = true:  Find nth prime >= nth_prime_approx.
/// @sieve_forward = false: Find nth prime <= nth_prime_approx.
///
template <bool sieve_forward, typename T>
T nth_prime_sieve2(uint64_t n,
                   T nth_prime_approx,
                   int max_threads)
{
  ASSERT(n > 0);
  ASSERT(nth_prime_approx > 0);

  // Calculate dist_approx such that:
  // nth_prime < nth_prime_approx + dist_approx
  double logx = std::log((double) nth_prime_approx);
  double log10_n = std::log10(max(n, 10));
  uint64_t dist_approx = (uint64_t) std::ceil(n * (logx + 1));
  dist_approx += uint64_t(logx * logx * log10_n);

  uint64_t root2 = (uint64_t) iroot<2>(nth_prime_approx);
  uint64_t root3 = (uint64_t) iroot<3>(nth_prime_approx);
  uint64_t max_thread_dist = uint64_t(root3 * 30);
  uint64_t thread_dist = in_between(240u, dist_approx, max_thread_dist);

  int main_threads = ideal_num_threads(dist_approx, max_threads, thread_dist);
  int max_threads_per_segment = in_between(1, ceil_div(max_threads, main_threads), 32);
  int threads_per_segment = ideal_num_threads(root2, max_threads_per_segment, thread_threshold);
  int total_threads = in_between(1, main_threads * threads_per_segment, max_threads);

  // Our nth_prime_sieve2 uses atomic memory accesses because
  // multiple threads process the same segment simultaneously.
  // When using a single thread per segment, we use the faster
  // nth_prime_sieve1 without atomic memory accesses.
  if (threads_per_segment == 1)
    return nth_prime_sieve1<sieve_forward>(n, nth_prime_approx, max_threads);

  aligned_vector<NthPrimeSieve2<T>> sieves(main_threads);

  T nth_prime = 0;
  uint64_t count = 0;
  uint64_t while_iters = 0;
  bool finished = false;
  double time;

  if (is_print())
  {
    print("");
    print("=== nth_prime_sieve ===");
    print_nth_prime_sieve(n, sieve_forward, nth_prime_approx, dist_approx,
        thread_dist, main_threads, threads_per_segment);
    time = get_time();
  }

  #pragma omp parallel num_threads(total_threads)
  #pragma omp single nowait
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

        #pragma omp task shared(sieves)
        {
          // Sieve the current segment [low, high].
          // If possible use fast 64-bit integer division
          // instead of slow 128-bit integer division.
          if ( low <= pstd::numeric_limits<uint64_t>::max() &&
              high <= pstd::numeric_limits<uint64_t>::max())
            sieves[t].sieve(uint64_t(low), uint64_t(high), threads_per_segment);
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

  if_unlikely(!nth_prime)
    throw primecount_error("Failed to find nth prime!");

  if (is_print())
    print_seconds(get_time() - time);

  return nth_prime;
}

#endif

template <typename T>
T nth_prime_sieve(T n,
                  T nth_prime_approx,
                  T count_approx,
                  int threads)
{
#if _OPENMP >= 201307

  // Use nth_prime_sieve2() if the CPU supports fast,
  // lock-free 64-bit atomics. Otherwise, fall back
  // to nth_prime_sieve1() to prevent the performance
  // penalties of software-emulated atomic accesses.
  #if __cplusplus >= 201703L
    constexpr bool atomic_uint64_lock_free = std::atomic<uint64_t>::is_always_lock_free;
  #else
    std::atomic<uint64_t> atomic_uint64(0);
    bool atomic_uint64_lock_free = atomic_uint64.is_lock_free();
  #endif

  if (atomic_uint64_lock_free && threads > 1)
  {
    if (count_approx < n)
      return nth_prime_sieve2<true>(uint64_t(n - count_approx),
                                    nth_prime_approx + 1,
                                    threads);
    else
      return nth_prime_sieve2<false>(uint64_t(1 + count_approx - n),
                                     nth_prime_approx,
                                     threads);
  }
#endif

  if (count_approx < n)
    return nth_prime_sieve1<true>(uint64_t(n - count_approx),
                                  nth_prime_approx + 1,
                                  threads);
  else
    return nth_prime_sieve1<false>(uint64_t(1 + count_approx - n),
                                   nth_prime_approx,
                                   threads);
}

} // namespace

namespace primecount {

int64_t nth_prime_sieve(int64_t n,
                        int64_t nth_prime_approx,
                        int64_t count_approx,
                        int threads)
{
  return ::nth_prime_sieve(n, nth_prime_approx, count_approx, threads);
}

#if defined(HAVE_INT128_T)

int128_t nth_prime_sieve(int128_t n,
                         int128_t nth_prime_approx,
                         int128_t count_approx,
                         int threads)
{
  return ::nth_prime_sieve(n, nth_prime_approx, count_approx, threads);
}

#endif

} // namespace
