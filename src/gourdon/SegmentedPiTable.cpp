///
/// @file  SegmentedPiTable.cpp
/// @brief The A and C formulas in Xavier Gourdon's prime counting
///        algorithm require looking up PrimePi[x] values with
///        x < x^(1/2). Since a PrimePi[x] lookup table of size x^(1/2)
///        would use too much memory we need a segmented PrimePi[x]
///        lookup table that uses only O(x^(1/4)) memory.
///
///        The SegmentedPiTable class is a compressed lookup table of
///        prime counts. Since the size of SegmentedPiTable is very
///        small and will always fit into the CPU's cache, we don't
///        use a bit array with maximum compression because this adds
///        significant overhead. Instead we use a bit array where each
///        bit corresponds to an odd integer. This compression scheme
///        provides very fast access since the bit array index can be
///        calculated using a single right shift instruction.
///
///        The algorithm of the easy special leaves and the usage of
///        the SegmentedPiTable are described in more detail in:
///        https://github.com/kimwalisch/primecount/blob/master/doc/Easy-Special-Leaves.pdf
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "SegmentedPiTable.hpp"

#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <imath.hpp>
#include <macros.hpp>
#include <min.hpp>

#include <stdint.h>
#include <algorithm>


namespace {

constexpr uint64_t bitmask(uint64_t n)
{
  return ((n + 1) / 2 == 64) ? 0xffffffffffffffffull
         : (1ull << ((n + 1) / 2)) - 1;
}

} // namespace

namespace primecount {

const Array<uint64_t, 128> SegmentedPiTable::unset_larger_ =
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

void SegmentedPiTable::init(uint64_t low, uint64_t high)
{
  ASSERT(low < high);
  ASSERT(low % 128 == 0);
  int threads = 1;
  uint64_t pi_low = 0;

  // In order to make the threads completely independent from
  // each other, each thread needs to compute PrimePi[low]
  // at the start of each newly assigned segment from the
  // LoadBalancer. However if a thread processes consecutive
  // segments, then we can compute PrimePi[low] in O(1) by
  // getting that value from the previous segment.
  if (low >= 2)
  {
    if (low == high_)
      pi_low = operator[](low - 1);
    else
      pi_low = pi_noprint(low - 1, threads);
  }

  low_ = low;
  high_ = high;
  uint64_t segment_size = high - low;
  uint64_t size = ceil_div(segment_size, 128);

  pi_.clear();
  pi_.resize(size);
  std::fill(pi_.begin(), pi_.end(), pi_t{0, 0});

  init_bits();
  init_count(pi_low);
}

/// Init pi[x] lookup table for [low, high[
void SegmentedPiTable::init_bits()
{
  // Iterate over primes >= 3
  uint64_t low = max(low_, 3);
  if (low >= high_)
    return;

  primesieve::iterator it(low, high_);
  uint64_t prime = 0;

  // For each prime in [low, high[ set the
  // corresponding bit in the pi[x] lookup table.
  while ((prime = it.next_prime()) < high_)
  {
    uint64_t i = (prime - low_) / 128;
    uint64_t prime_bit = 1ull << (prime % 128 / 2);
    pi_[i].bits |= prime_bit;
  }
}

void SegmentedPiTable::init_count(uint64_t pi_low)
{
  // Workaround needed for prime 2 since
  // we are sieving with primes >= 3.
  if (low_ < 2 && high_ > 2)
    pi_low += 1;

  uint64_t j = ceil_div((high_ - low_), 128);

  // Count 1 bits (primes) in pi[x] lookup table
  for (uint64_t i = 0; i < j; i++)
  {
    pi_[i].count = pi_low;
    pi_low += popcnt64(pi_[i].bits);
  }
}

} // namespace
