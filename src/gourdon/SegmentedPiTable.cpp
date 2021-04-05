///
/// @file  SegmentedPiTable.cpp
/// @brief The A and C formulas in Xavier Gourdon's prime counting
///        algorithm require looking up PrimePi[n] values with
///        n < x^(1/2). Since a PrimePi[n] lookup table of size x^(1/2)
///        would use too much memory we need a segmented PrimePi[n]
///        lookup table that uses only O(x^(1/3)) memory.
///
///        The SegmentedPiTable class is a compressed lookup table of
///        prime counts. Each bit of the lookup table corresponds to
///        an integer that is not divisible by 2, 3 and 5. The 8 bits
///        of each byte correspond to the offsets { 1, 7, 11, 13, 17,
///        19, 23, 29 }. Since our lookup table uses the uint64_t data
///        type, one array element (8 bytes) corresponds to an
///        interval of size 30 * 8 = 240.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <SegmentedPiTable.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <imath.hpp>
#include <min.hpp>

#include <stdint.h>
#include <cstring>

namespace primecount {

SegmentedPiTable::SegmentedPiTable(uint64_t max_high,
                                   uint64_t segment_size,
                                   int threads)
  : counts_(threads),
    max_high_(max_high),
    threads_(threads)
{
  // The threads in our AC algorithm are not completely
  // independent from each other. After each segment all
  // threads need to be synchronized. On servers with a
  // large number of CPU cores this can add a lot of
  // overhead. For this reason we set a large minimum
  // segment size here (2 MiB, should be >= CPU L2 cache
  // and <= CPU L3 cache) to avoid such scaling issues.
  uint64_t numbers_per_byte = 240 / sizeof(pi_t);
  uint64_t min_segment_size = (2 << 20) * numbers_per_byte;
  segment_size_ = max(segment_size, min_segment_size);
  segment_size_ = min(segment_size_, max_high_);

  // In order to simplify multi-threading we set low,
  // high and segment_size % 240 == 0.
  segment_size_ += 240 - segment_size_ % 240;
  pi_.resize(segment_size_ / 240);
}

/// Iterate over the primes inside the segment [low, high[
/// and initialize the pi[x] lookup table. The pi[x]
/// lookup table returns the number of primes <= x for
/// low <= x < high.
///
void SegmentedPiTable::init(uint64_t low, uint64_t high)
{
  #pragma omp master
  {
    if (low > 0)
      pi_low_ = operator[](low - 1);
    low_ = low;
    high_ = min(high, max_high_);
  }

  uint64_t thread_size = segment_size_ / threads_;
  uint64_t min_thread_size = (uint64_t) 1e6;
  thread_size = max(min_thread_size, thread_size);
  thread_size += 240 - thread_size % 240;

  #pragma omp for
  for (int t = 0; t < threads_; t++)
  {
    uint64_t thread_low = low + thread_size * t;
    uint64_t thread_high = thread_low + thread_size;
    thread_high = min(thread_high, high);

    if (thread_low < thread_high)
      init_bits(low, thread_low, thread_high, t);
  }

  #pragma omp for
  for (int t = 0; t < threads_; t++)
  {
    uint64_t thread_low = low + thread_size * t;
    uint64_t thread_high = thread_low + thread_size;
    thread_high = min(thread_high, high);

    if (thread_low < thread_high)
      init_count(low, thread_low, thread_high, t);
  }
}

/// Each thread computes PrimePi [thread_low, thread_high[
void SegmentedPiTable::init_bits(uint64_t low,
                                 uint64_t thread_low,
                                 uint64_t thread_high,
                                 uint64_t thread_num)
{
  // Zero initialize pi vector
  uint64_t i = (thread_low - low) / 240;
  uint64_t j = ceil_div((thread_high - low), 240);
  std::memset(&pi_[i], 0, (j - i) * sizeof(pi_t));

  // Iterate over primes > 5
  primesieve::iterator it(max(thread_low, 5), thread_high);
  uint64_t count = 0;
  uint64_t prime = 0;

  // Each thread iterates over the primes
  // inside [thread_low, thread_high[ and initializes
  // the pi[x] lookup table.
  while ((prime = it.next_prime()) < thread_high)
  {
    uint64_t p = prime - low;
    pi_[p / 240].bits |= set_bit_[p % 240];
    count += 1;
  }

  counts_[thread_num] = count;
}

/// Each thread computes PrimePi [thread_low, thread_high[
void SegmentedPiTable::init_count(uint64_t low,
                                  uint64_t thread_low,
                                  uint64_t thread_high,
                                  uint64_t thread_num)
{
  // First compute PrimePi[thread_low - 1]
  uint64_t count = pi_low_;
  for (uint64_t i = 0; i < thread_num; i++)
    count += counts_[i];

  // Convert to array indexes
  uint64_t i = (thread_low - low) / 240;
  uint64_t stop_idx = ceil_div((thread_high - low), 240);

  for (; i < stop_idx; i++)
  {
    pi_[i].count = count;
    count += popcnt64(pi_[i].bits);
  }
}

} // namespace
