///
/// @file  BitSieve128.cpp
/// @brief The BitSieve128 base class contains lookup tables that are
///        needed to implement a prime sieving algorithm where each bit
///        corresponds to an odd integer. The BitSieve128 class uses
///        the uint64_t data type for its sieve array, hence each sieve
///        array element corresponds to an interval of 64 * 2 = 128.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <BitSieve128.hpp>

#include <stdint.h>
#include <array>

namespace {

constexpr uint64_t bitmask(uint64_t n)
{
  return ((n + 1) / 2 == 64) ? 0xffffffffffffffffull
         : (1ull << ((n + 1) / 2)) - 1;
}

} // namespace

namespace primecount {

/// Bitmasks needed to unset a specific bit in the sieve array
const std::array<uint64_t, 128> BitSieve128::unset_bit_ =
{
  ~(1ull << 0), ~(1ull << 0),
  ~(1ull << 1), ~(1ull << 1),
  ~(1ull << 2), ~(1ull << 2),
  ~(1ull << 3), ~(1ull << 3),
  ~(1ull << 4), ~(1ull << 4),
  ~(1ull << 5), ~(1ull << 5),
  ~(1ull << 6), ~(1ull << 6),
  ~(1ull << 7), ~(1ull << 7),
  ~(1ull << 8), ~(1ull << 8),
  ~(1ull << 9), ~(1ull << 9),
  ~(1ull << 10), ~(1ull << 10),
  ~(1ull << 11), ~(1ull << 11),
  ~(1ull << 12), ~(1ull << 12),
  ~(1ull << 13), ~(1ull << 13),
  ~(1ull << 14), ~(1ull << 14),
  ~(1ull << 15), ~(1ull << 15),
  ~(1ull << 16), ~(1ull << 16),
  ~(1ull << 17), ~(1ull << 17),
  ~(1ull << 18), ~(1ull << 18),
  ~(1ull << 19), ~(1ull << 19),
  ~(1ull << 20), ~(1ull << 20),
  ~(1ull << 21), ~(1ull << 21),
  ~(1ull << 22), ~(1ull << 22),
  ~(1ull << 23), ~(1ull << 23),
  ~(1ull << 24), ~(1ull << 24),
  ~(1ull << 25), ~(1ull << 25),
  ~(1ull << 26), ~(1ull << 26),
  ~(1ull << 27), ~(1ull << 27),
  ~(1ull << 28), ~(1ull << 28),
  ~(1ull << 29), ~(1ull << 29),
  ~(1ull << 30), ~(1ull << 30),
  ~(1ull << 31), ~(1ull << 31),
  ~(1ull << 32), ~(1ull << 32),
  ~(1ull << 33), ~(1ull << 33),
  ~(1ull << 34), ~(1ull << 34),
  ~(1ull << 35), ~(1ull << 35),
  ~(1ull << 36), ~(1ull << 36),
  ~(1ull << 37), ~(1ull << 37),
  ~(1ull << 38), ~(1ull << 38),
  ~(1ull << 39), ~(1ull << 39),
  ~(1ull << 40), ~(1ull << 40),
  ~(1ull << 41), ~(1ull << 41),
  ~(1ull << 42), ~(1ull << 42),
  ~(1ull << 43), ~(1ull << 43),
  ~(1ull << 44), ~(1ull << 44),
  ~(1ull << 45), ~(1ull << 45),
  ~(1ull << 46), ~(1ull << 46),
  ~(1ull << 47), ~(1ull << 47),
  ~(1ull << 48), ~(1ull << 48),
  ~(1ull << 49), ~(1ull << 49),
  ~(1ull << 50), ~(1ull << 50),
  ~(1ull << 51), ~(1ull << 51),
  ~(1ull << 52), ~(1ull << 52),
  ~(1ull << 53), ~(1ull << 53),
  ~(1ull << 54), ~(1ull << 54),
  ~(1ull << 55), ~(1ull << 55),
  ~(1ull << 56), ~(1ull << 56),
  ~(1ull << 57), ~(1ull << 57),
  ~(1ull << 58), ~(1ull << 58),
  ~(1ull << 59), ~(1ull << 59),
  ~(1ull << 60), ~(1ull << 60),
  ~(1ull << 61), ~(1ull << 61),
  ~(1ull << 62), ~(1ull << 62),
  ~(1ull << 63), ~(1ull << 63)
};

/// unset_larger_[x % 128] returns a bitmask where the bits
/// corresponding to numbers > x % 128 have been turned off
/// (while all the other bits are left unchanged).
/// 
const std::array<uint64_t, 128> BitSieve128::unset_larger_ =
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

} // namespace
