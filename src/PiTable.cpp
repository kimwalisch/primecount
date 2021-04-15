///
/// @file  PiTable.cpp
/// @brief The PiTable class is a compressed lookup table of prime
///        counts. Each bit of the lookup table corresponds to an
///        integer that is not divisible by 2, 3 and 5. The 8 bits of
///        each byte correspond to the offsets { 1, 7, 11, 13, 17, 19,
///        23, 29 }. Since our lookup table uses the uint64_t data
///        type, one array element (8 bytes) corresponds to an
///        interval of size 30 * 8 = 240.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <imath.hpp>
#include <min.hpp>

#include <stdint.h>
#include <cstring>

namespace primecount {

PiTable::PiTable(uint64_t limit, int threads) :
  limit_(limit)
{
  uint64_t size = limit + 1;
  uint64_t thread_threshold = (uint64_t) 1e7;
  threads = ideal_num_threads(threads, size, thread_threshold);
  uint64_t thread_size = size / threads;
  thread_size = max(thread_threshold, thread_size);
  thread_size += 240 - thread_size % 240;
  pi_.resize(ceil_div(size, 240));
  counts_.resize(threads);

  #pragma omp parallel num_threads(threads)
  {
    #pragma omp for
    for (int t = 0; t < threads; t++)
    {
      uint64_t start = thread_size * t;
      uint64_t stop = start + thread_size;
      stop = min(stop, size);

      if (start < stop)
        init_bits(start, stop, t);
    }

    #pragma omp for
    for (int t = 0; t < threads; t++)
    {
      uint64_t start = thread_size * t;
      uint64_t stop = start + thread_size;
      stop = min(stop, size);

      if (start < stop)
        init_count(start, stop, t);
    }
  }
}

/// Each thread computes PrimePi [start, stop[
void PiTable::init_bits(uint64_t start,
                        uint64_t stop,
                        uint64_t thread_num)
{
  // Zero initialize pi vector
  uint64_t i = start / 240;
  uint64_t j = ceil_div(stop, 240);
  std::memset(&pi_[i], 0, (j - i) * sizeof(pi_t));

  // Iterate over primes > 5
  start = max(start, 5);
  primesieve::iterator it(start, stop);
  uint64_t count = 0;
  uint64_t prime = 0;

  while ((prime = it.next_prime()) < stop)
  {
    uint64_t prime_bit = set_bit_[prime % 240];
    pi_[prime / 240].bits |= prime_bit;
    count += 1;
  }

  counts_[thread_num] = count;
}

/// Each thread computes PrimePi [start, stop[
void PiTable::init_count(uint64_t start,
                         uint64_t stop,
                         uint64_t thread_num)
{
  // First compute PrimePi[start - 1]
  uint64_t count = pi_tiny_[5];
  for (uint64_t i = 0; i < thread_num; i++)
    count += counts_[i];

  // Convert to array indexes
  uint64_t i = start / 240;
  uint64_t stop_idx = ceil_div(stop, 240);

  for (; i < stop_idx; i++)
  {
    pi_[i].count = count;
    count += popcnt64(pi_[i].bits);
  }
}

} // namespace
