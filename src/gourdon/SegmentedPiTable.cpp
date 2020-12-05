///
/// @file  SegmentedPiTable.cpp
/// @brief The A and C formulas in Xavier Gourdon's prime counting
///        algorithm require looking up PrimePi[n] values with
///        n < x^(1/2). Since a PrimePi[n] lookup table of size
///        x^(1/2) would use too much memory we need a segmented
///        PrimePi[n] lookup table that uses only O(z) memory.
///
///        The SegmentedPiTable is based on the PiTable class which
///        is a compressed lookup table for prime counts. Each bit
///        in the lookup table corresponds to an odd integer and that
///        bit is set to 1 if the integer is a prime. PiTable uses
///        only (n / 8) bytes of memory and returns the number of
///        primes <= n in O(1) operations.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
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
#include <array>
#include <cstring>

namespace {

constexpr uint64_t bitmask(uint64_t n)
{
  return ((n + 1) / 2 == 64) ? 0xffffffffffffffffull
         : (1ull << ((n + 1) / 2)) - 1;
}

} // namespace

namespace primecount {

const std::array<uint64_t, 128> SegmentedPiTable::unset_bits_ =
{
  bitmask(0),   bitmask(1),   bitmask(2),   bitmask(3),
  bitmask(4),   bitmask(5),   bitmask(6),   bitmask(7),
  bitmask(8),   bitmask(9),   bitmask(10),  bitmask(11),
  bitmask(12),  bitmask(13),  bitmask(14),  bitmask(15),
  bitmask(16),  bitmask(17),  bitmask(18),  bitmask(19),
  bitmask(20),  bitmask(21),  bitmask(22),  bitmask(23),
  bitmask(24),  bitmask(25),  bitmask(26),  bitmask(27),
  bitmask(28),  bitmask(29),  bitmask(30),  bitmask(31),
  bitmask(32),  bitmask(33),  bitmask(34),  bitmask(35),
  bitmask(36),  bitmask(37),  bitmask(38),  bitmask(39),
  bitmask(40),  bitmask(41),  bitmask(42),  bitmask(43),
  bitmask(44),  bitmask(45),  bitmask(46),  bitmask(47),
  bitmask(48),  bitmask(49),  bitmask(50),  bitmask(51),
  bitmask(52),  bitmask(53),  bitmask(54),  bitmask(55),
  bitmask(56),  bitmask(57),  bitmask(58),  bitmask(59),
  bitmask(60),  bitmask(61),  bitmask(62),  bitmask(63),
  bitmask(64),  bitmask(65),  bitmask(66),  bitmask(67),
  bitmask(68),  bitmask(69),  bitmask(70),  bitmask(71),
  bitmask(72),  bitmask(73),  bitmask(74),  bitmask(75),
  bitmask(76),  bitmask(77),  bitmask(78),  bitmask(79),
  bitmask(80),  bitmask(81),  bitmask(82),  bitmask(83),
  bitmask(84),  bitmask(85),  bitmask(86),  bitmask(87),
  bitmask(88),  bitmask(89),  bitmask(90),  bitmask(91),
  bitmask(92),  bitmask(93),  bitmask(94),  bitmask(95),
  bitmask(96),  bitmask(97),  bitmask(98),  bitmask(99),
  bitmask(100), bitmask(101), bitmask(102), bitmask(103),
  bitmask(104), bitmask(105), bitmask(106), bitmask(107),
  bitmask(108), bitmask(109), bitmask(110), bitmask(111),
  bitmask(112), bitmask(113), bitmask(114), bitmask(115),
  bitmask(116), bitmask(117), bitmask(118), bitmask(119),
  bitmask(120), bitmask(121), bitmask(122), bitmask(123),
  bitmask(124), bitmask(125), bitmask(126), bitmask(127)
};

SegmentedPiTable::SegmentedPiTable(uint64_t low,
                                   uint64_t limit,
                                   uint64_t segment_size,
                                   int threads)
  : counts_(threads),
    low_(low),
    max_high_(limit + 1),
    threads_(threads)
{
  // Each bit of the pi[x] lookup table corresponds
  // to an odd integer, so there are 16 numbers per
  // byte. However we also store 64-bit prime_count
  // values in the pi[x] lookup table, hence each byte
  // only corresponds to 8 numbers.
  uint64_t numbers_per_byte = 8;

  // Minimum segment size = 256 KiB (L2 cache size),
  // a large segment size improves load balancing.
  uint64_t min_segment_size = 256 * (1 << 10) * numbers_per_byte;
  segment_size_ = max(segment_size, min_segment_size);
  segment_size_ = min(segment_size_, max_high_);

  // In order to simplify multi-threading we set low,
  // high and segment_size % 128 == 0.
  segment_size_ += 128 - segment_size_ % 128;
  int thread_threshold = (int) 1e7;
  threads_ = ideal_num_threads(threads, segment_size_, thread_threshold);

  high_ = low_ + segment_size_;
  high_ = std::min(high_, max_high_);
  pi_.resize(segment_size_ / 128);

  low = max(low, 1);
  pi_low_ = pi_simple(low - 1, threads);
  init();
}

/// Increase low & high and initialize the next segment.
void SegmentedPiTable::next()
{
  // pi_low_ must be initialized before updating the
  // member variables for the next segment.
  pi_low_ = operator[](high_ - 1);

  low_ = high_;
  high_ = low_ + segment_size_;
  high_ = std::min(high_, max_high_);

  if (finished())
    return;

  init();
}

/// Iterate over the primes inside the segment [low, high[
/// and initialize the pi[x] lookup table. The pi[x]
/// lookup table returns the number of primes <= x for
/// low <= x < high.
///
void SegmentedPiTable::init()
{
  uint64_t thread_size = segment_size_ / threads_;
  uint64_t min_thread_size = (uint64_t) 1e7;
  thread_size = max(min_thread_size, thread_size);
  thread_size += 128 - thread_size % 128;

  #pragma omp parallel num_threads(threads_)
  {
    #pragma omp for schedule(dynamic)
    for (int t = 0; t < threads_; t++)
    {
      uint64_t start = low_ + thread_size * t;
      uint64_t stop = start + thread_size;
      stop = min(stop, high_);

      if (start < stop)
        init_bits(start, stop, t);
    }

    #pragma omp for schedule(dynamic)
    for (int t = 0; t < threads_; t++)
    {
      uint64_t start = low_ + thread_size * t;
      uint64_t stop = start + thread_size;
      stop = min(stop, high_);

      if (start < stop)
        init_prime_count(start, stop, t);
    }
  }
}

/// Each thread computes PrimePi [start, stop[
void SegmentedPiTable::init_bits(uint64_t start,
                                 uint64_t stop,
                                 uint64_t thread_num)
{
  // Zero initialize pi vector
  uint64_t i = (start - low_) / 128;
  uint64_t j = ceil_div((stop - low_), 128);
  std::memset(&pi_[i], 0, (j - i) * sizeof(pi_t));

  // Since we store only odd numbers in our lookup table,
  // we cannot store 2 which is the only even prime.
  // As a workaround we mark 1 as a prime (1st bit) and
  // add a check to return 0 for pi[1].
  if (start <= 1)
    pi_[0].bits |= 1;

  // Iterate over primes > 2
  primesieve::iterator it(max(start, 2), stop);
  uint64_t count = (start <= 2);
  uint64_t prime = 0;

  // Each thread iterates over the primes
  // inside [start, stop[ and initializes
  // the pi[x] lookup table.
  while ((prime = it.next_prime()) < stop)
  {
    uint64_t p = prime - low_;
    uint64_t prime_bit = 1ull << (p % 128 / 2);
    pi_[p / 128].bits |= prime_bit;
    count += 1;
  }

  counts_[thread_num] = count;
}

/// Each thread computes PrimePi [start, stop[
void SegmentedPiTable::init_prime_count(uint64_t start,
                                        uint64_t stop,
                                        uint64_t thread_num)
{
  // First compute PrimePi[start - 1]
  uint64_t count = pi_low_;
  for (uint64_t i = 0; i < thread_num; i++)
    count += counts_[i];

  // Convert to array indexes
  uint64_t i = (start - low_) / 128;
  uint64_t stop_idx = ceil_div((stop - low_), 128);

  for (; i < stop_idx; i++)
  {
    pi_[i].prime_count = count;
    count += popcnt64(pi_[i].bits);
  }
}

} // namespace
