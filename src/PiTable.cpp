///
/// @file  PiTable.cpp
/// @brief The PiTable class is a compressed lookup table for
///        prime counts. Each bit in the lookup table corresponds
///        to an odd integer and that bit is set to 1 if the
///        integer is a prime. PiTable uses only (n / 8) bytes of
///        memory and returns the number of primes <= n in O(1)
///        operations.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <aligned_vector.hpp>
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
  thread_size += 128 - thread_size % 128;
  pi_.resize(ceil_div(size, 128));
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
        init_prime_count(start, stop, t);
    }
  }
}

/// Each thread computes PrimePi [start, stop[
void PiTable::init_bits(uint64_t start,
                        uint64_t stop,
                        uint64_t thread_num)
{
  // Zero initialize pi vector
  uint64_t i = start / 128;
  uint64_t j = ceil_div(stop, 128);
  std::memset(&pi_[i], 0, (j - i) * sizeof(pi_t));

  // Iterate over primes > 2
  primesieve::iterator it(max(start, 2), stop);
  uint64_t count = 0;
  uint64_t prime = 0;

  while ((prime = it.next_prime()) < stop)
  {
    uint64_t prime_bit = set_bit_[prime % 128];
    pi_[prime / 128].bits |= prime_bit;
    count += 1;
  }

  counts_[thread_num] = count;
}

/// Each thread computes PrimePi [start, stop[
void PiTable::init_prime_count(uint64_t start,
                               uint64_t stop,
                               uint64_t thread_num)
{
  // First compute PrimePi[start - 1]
  uint64_t count = pi_tiny_[2];
  for (uint64_t i = 0; i < thread_num; i++)
    count += counts_[i];

  // Convert to array indexes
  uint64_t i = start / 128;
  uint64_t stop_idx = ceil_div(stop, 128);

  for (; i < stop_idx; i++)
  {
    pi_[i].prime_count = count;
    count += popcnt64(pi_[i].bits);
  }
}

} // namespace
