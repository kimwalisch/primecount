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

constexpr uint64_t min_sieving_prime_dist = UINT64_C(100000000);
constexpr int max_thread_group_size = 8;

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

  template <typename UT>
  void sieve_parallel(UT low, UT high, int max_threads)
  {
    init(low, high);
    int active_threads = get_active_threads(max_threads);

    // Let the current segment task process one partition itself
    // and spawn only the remaining helper tasks.
    #pragma omp taskgroup
    {
      for (int thread_id = 1; thread_id < active_threads; thread_id++)
      {
        #pragma omp task firstprivate(thread_id, active_threads)
        cross_off_prime_range(thread_id, active_threads);
      }

      cross_off_prime_range(0, active_threads);
    }

    count_primes();
  }

  void cross_off_prime_range(int thread_id, int threads)
  {
    using UT = typename pstd::make_unsigned<T>::type;

    if (!size_)
      return;

    uint64_t prime_start = sqrt_high_ * thread_id / threads;
    uint64_t prime_stop = sqrt_high_ * (thread_id + 1) / threads;
    prime_start += (thread_id > 0);
    prime_start = max(prime_start, UINT64_C(7));

    if (prime_start > prime_stop)
      return;

    primesieve::iterator iter(prime_start, prime_stop);
    uint64_t prime;
    UT low = (UT) low_;
    auto* sieve = sieve_.get();

    if (use_prime_squares_)
    {
      while ((prime = iter.next_prime()) <= prime_stop)
      {
        // Calculate first multiple > low
        // that is not divisible by 2.
        UT q = fast_div(low, prime);
        UT n = prime * (q + 1 + (q & 1));
        n = max(n, UT(prime) * prime);
        ASSERT(n > prime);
        ASSERT(n % 2 != 0);
        uint64_t i = (uint64_t) (n - low);

        // Cross-off multiples
        for (; i <= i_max_; i += prime * 2)
          sieve[i / 240].fetch_and(unset_bit_[i % 240], std::memory_order_relaxed);
      }
    }
    else
    {
      while ((prime = iter.next_prime()) <= prime_stop)
      {
        // Calculate first multiple > low
        // that is not divisible by 2.
        UT q = fast_div(low, prime);
        UT n = prime * (q + 1 + (q & 1));
        ASSERT(n > prime);
        ASSERT(n % 2 != 0);
        uint64_t i = (uint64_t) (n - low);

        // Cross-off multiples
        for (; i <= i_max_; i += prime * 2)
          sieve[i / 240].fetch_and(unset_bit_[i % 240], std::memory_order_relaxed);
      }
    }
  }

  int get_active_threads(int max_threads) const
  {
    int threads = (int) (sqrt_high_ / min_sieving_prime_dist);
    return in_between(1, threads, max_threads);
  }

  void count_primes()
  {
    uint64_t count = 0;

    // Count primes (1 bits)
    for (uint64_t i = 0; i < size_; i++)
      count += popcnt64(load_word(i));

    count_ = count;
  }

  T find_nth_prime(uint64_t n) const
  {
    ASSERT(n > 0);
    ASSERT(n <= count_);

    uint64_t count = 0;

    for (uint64_t i = 0; i < size_; i++)
    {
      uint64_t bits = load_word(i);
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
  template <typename UT>
  void init(UT low, UT high)
  {
    if (high < 2)
    {
      low_ = 0;
      count_ = 0;
      size_ = 0;
      sqrt_high_ = 0;
      i_max_ = 0;
      use_prime_squares_ = false;
      return;
    }

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
    i_max_ = (uint64_t) (high - low);
    size_ = (uint64_t) ceil_div(i_max_ + 1, 240);
    sqrt_high_ = (uint64_t) isqrt(high);
    use_prime_squares_ = sqrt_high_ >= low;
    count_ = 0;
    uint64_t first_word = ~0ull;
    uint64_t last_word = ~0ull;
    first_word &= unset_smaller_[old_low % 240];
    last_word &= unset_larger_[high % 240];

    if (size_ > sieve_size_)
    {
      sieve_.reset(new std::atomic<uint64_t>[size_]);
      sieve_size_ = size_;
    }

    auto* sieve = sieve_.get();
    uint64_t word = first_word;

    if (size_ == 1)
      word &= last_word;

    sieve[0].store(word, std::memory_order_relaxed);

    for (uint64_t i = 1; i + 1 < size_; i++)
      sieve[i].store(~0ull, std::memory_order_relaxed);

    if (size_ > 1)
      sieve[size_ - 1].store(~0ull & last_word, std::memory_order_relaxed);
  }

  ALWAYS_INLINE uint64_t load_word(uint64_t i) const
  {
    return sieve_[i].load(std::memory_order_relaxed);
  }

  T low_ = 0;
  uint64_t count_ = 0;
  uint64_t size_ = 0;
  uint64_t sqrt_high_ = 0;
  uint64_t i_max_ = 0;
  bool use_prime_squares_ = false;
  std::unique_ptr<std::atomic<uint64_t>[]> sieve_;
  uint64_t sieve_size_ = 0;
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
  int max_threads = threads;

  threads = ideal_num_threads(dist_approx, threads, thread_dist);
  int thread_group_size = 1;

#ifdef _OPENMP
  {
    using UT = typename pstd::make_unsigned<T>::type;
    uint64_t sqrt_n = (uint64_t) isqrt((UT) nth_prime_approx);
    int thread_budget = max(1, max_threads / threads);
    thread_group_size = (int) (sqrt_n / min_sieving_prime_dist);
    thread_group_size = in_between(1, thread_group_size, max_thread_group_size);
    thread_group_size = std::min(thread_group_size, thread_budget);
  }
#endif

  int total_threads = threads * thread_group_size;
  aligned_vector<NthPrimeSieve<T>> sieves(threads);
  bool print_vars = is_print();
  bool finished = false;
  double time;

  if (print_vars)
  {
    print("");
    print("=== nth_prime_sieve ===");
    print_nth_prime_sieve(n, sieve_forward, nth_prime_approx, dist_approx, thread_dist, total_threads);
    if (thread_group_size > 1)
    {
      print("segment_threads", threads);
      print("thread_group_size", thread_group_size);
    }
    time = get_time();
  }

  #pragma omp parallel num_threads(total_threads)
  {
    #pragma omp single nowait
    while (!finished)
    {
      uint64_t current_iter = while_iters++;

      #pragma omp taskgroup
      {
        for (int t = 0; t < threads; t++)
        {
          // Unsigned integer division is usually
          // faster than signed integer division.
          using UT = typename pstd::make_unsigned<T>::type;
          uint64_t i = current_iter * threads + t;
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

          #pragma omp task firstprivate(t, low, high, thread_group_size) shared(sieves)
          {
            // Sieve the current segment [low, high].
            // If possible use fast 64-bit bit integer division
            // instead of slow 128-bit integer division.
            if ( low <= pstd::numeric_limits<uint64_t>::max() &&
                high <= pstd::numeric_limits<uint64_t>::max())
              sieves[t].sieve_parallel((uint64_t) low, (uint64_t) high, thread_group_size);
            else
              sieves[t].sieve_parallel(low, high, thread_group_size);
          }
        }
      }

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
