///
/// @file  SegmentedPiTable.cpp
/// @brief The A and C formulas in Xavier Gourdon's prime counting
///        algorithm require looking up PrimePi[n] values with
///        n <= x^(1/2). Since a PrimePi[n] lookup table of size
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
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <SegmentedPiTable.hpp>
#include <primesieve.hpp>

#include <stdint.h>
#include <array>
#include <vector>

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

SegmentedPiTable::SegmentedPiTable(uint64_t sqrtx,
                                   uint64_t segment_size)
  : max_(sqrtx + 1)
{
  // Each bit of the pi[x] lookup table corresponds
  // to an odd integer, so there are 16 numbers per
  // byte. However we also store 64-bit prime_count
  // values in the pi[x] lookup table, hence each byte
  // only corresponds to 8 numbers.
  uint64_t numbers_per_byte = 8;

  // Minimum segment size = 256 KiB (L2 cache size),
  // A large segment size improves load balancing.
  uint64_t min_segment_size = 256 * (1 << 10) * numbers_per_byte;
  segment_size_ = std::max(segment_size, min_segment_size);
  segment_size_ = std::min(segment_size_, max_);
  segment_size_ = segment_size + segment_size % 2;

  high_ = segment_size_;
  high_ = std::min(high_, max_);
  pi_.resize(segment_size_ / 128 + 1);
  primesieve::iterator it(2, high_);

  uint64_t pix = 0;
  uint64_t prime = 0;

  // Since we store only odd numbers in our lookup table,
  // we cannot store 2 which is the only even prime.
  // As a workaround we mark 1 as a prime (1st bit) and
  // add a check to return 0 for pi[1].
  if (low_ <= 1)
    pi_[0].bits = 1;

  while ((prime = it.next_prime()) < high_)
  {
    uint64_t p = prime - low_;
    pi_[p / 128].bits |= 1ull << (p % 128 / 2);
  }

  for (auto& i : pi_)
  {
    i.prime_count = pix;
    pix += popcnt64(i.bits);
  }
}

/// Initialize next segment to [high, high + segment_size[
/// (the current segment is [low, high[).
///
void SegmentedPiTable::next()
{
  assert(low_ <= max_);
  uint64_t prime = 0;
  uint64_t pix = operator[](high_ - 1);

  low_ = high_;
  high_ = low_ + segment_size_;
  high_ = std::min(high_, max_);

  if (finished())
    return;

  // Reset pi[x] lookup table
  for (auto& i : pi_)
  {
    i.bits = 0;
    i.prime_count = 0;
  }

  primesieve::iterator it(low_ - 1, high_);

  // Mark prime numbers
  while ((prime = it.next_prime()) < high_)
  {
    uint64_t p = prime - low_;
    pi_[p / 128].bits |= 1ull << (p % 128 / 2);
  }

  // Update prime counts
  for (auto& i : pi_)
  {
    i.prime_count = pix;
    pix += popcnt64(i.bits);
  }
}

} // namespace
