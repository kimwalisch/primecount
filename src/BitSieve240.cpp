///
/// @file  BitSieve240.cpp
/// @brief The BitSieve240 base class contains lookup tables that are
///        needed to implement a prime sieving algorithm where each
///        bit corresponds to an integer that is not divisible by 2,
///        3 and 5. The 8 bits of each byte correspond to the offsets
///        { 1, 7, 11, 13, 17, 19, 23, 29 }. Since the sieve array
///        uses the uint64_t data type, one sieve array element
///        (8 bytes) corresponds to an interval of size 30 * 8 = 240.
///
/// Copyright (C) 2023 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <BitSieve240.hpp>
#include <Vector.hpp>

#include <stdint.h>

namespace {

/// The 8 bits in each byte of the sieve array correspond
/// to the offsets { 1, 7, 11, 13, 17, 19, 23, 29 }.
///
constexpr int left_shift(int n)
{
  return   (n % 30 <= 1)  ? (n / 30 * 8) + 0
         : (n % 30 <= 7)  ? (n / 30 * 8) + 1
         : (n % 30 <= 11) ? (n / 30 * 8) + 2
         : (n % 30 <= 13) ? (n / 30 * 8) + 3
         : (n % 30 <= 17) ? (n / 30 * 8) + 4
         : (n % 30 <= 19) ? (n / 30 * 8) + 5
         : (n % 30 <= 23) ? (n / 30 * 8) + 6
         : (n / 30 * 8) + 7;
}

/// The 8 bits in each byte of the sieve array correspond
/// to the offsets { 1, 7, 11, 13, 17, 19, 23, 29 }.
///
constexpr int right_shift(int n)
{
  return   (n % 30 >= 29) ? 56 - (n / 30 * 8)
         : (n % 30 >= 23) ? 57 - (n / 30 * 8)
         : (n % 30 >= 19) ? 58 - (n / 30 * 8)
         : (n % 30 >= 17) ? 59 - (n / 30 * 8)
         : (n % 30 >= 13) ? 60 - (n / 30 * 8)
         : (n % 30 >= 11) ? 61 - (n / 30 * 8)
         : (n % 30 >= 7)  ? 62 - (n / 30 * 8)
         : (n % 30 >= 1)  ? 63 - (n / 30 * 8)
         : 64 - (n / 30 * 8);
}

/// Bitmask to unset bits <= n
constexpr uint64_t unset_s(int n)
{
  return ~0ull << left_shift(n);
}

/// Bitmask to unset bits >= n
constexpr uint64_t unset_l(int n)
{
  return (n == 0) ? 0 : ~0ull >> right_shift(n);
}

} // namespace

namespace primecount {

/// pi(x) for x < 6
const Array<uint64_t, 6> BitSieve240::pi_tiny_ = { 0, 0, 1, 2, 2, 3 };

const Array<uint64_t, 64> BitSieve240::bit_values_ =
{
  0 * 30 + 1, 0 * 30 + 7, 0 * 30 + 11, 0 * 30 + 13, 0 * 30 + 17, 0 * 30 + 19, 0 * 30 + 23, 0 * 30 + 29,
  1 * 30 + 1, 1 * 30 + 7, 1 * 30 + 11, 1 * 30 + 13, 1 * 30 + 17, 1 * 30 + 19, 1 * 30 + 23, 1 * 30 + 29,
  2 * 30 + 1, 2 * 30 + 7, 2 * 30 + 11, 2 * 30 + 13, 2 * 30 + 17, 2 * 30 + 19, 2 * 30 + 23, 2 * 30 + 29,
  3 * 30 + 1, 3 * 30 + 7, 3 * 30 + 11, 3 * 30 + 13, 3 * 30 + 17, 3 * 30 + 19, 3 * 30 + 23, 3 * 30 + 29,
  4 * 30 + 1, 4 * 30 + 7, 4 * 30 + 11, 4 * 30 + 13, 4 * 30 + 17, 4 * 30 + 19, 4 * 30 + 23, 4 * 30 + 29,
  5 * 30 + 1, 5 * 30 + 7, 5 * 30 + 11, 5 * 30 + 13, 5 * 30 + 17, 5 * 30 + 19, 5 * 30 + 23, 5 * 30 + 29,
  6 * 30 + 1, 6 * 30 + 7, 6 * 30 + 11, 6 * 30 + 13, 6 * 30 + 17, 6 * 30 + 19, 6 * 30 + 23, 6 * 30 + 29,
  7 * 30 + 1, 7 * 30 + 7, 7 * 30 + 11, 7 * 30 + 13, 7 * 30 + 17, 7 * 30 + 19, 7 * 30 + 23, 7 * 30 + 29
};

/// Bitmasks needed to set a specific bit in the sieve array
const Array<uint64_t, 240> BitSieve240::set_bit_ =
{
  0ull, 1ull << 0, 0ull, 0ull, 0ull,
  0ull, 0ull, 1ull << 1, 0ull, 0ull,
  0ull, 1ull << 2, 0ull, 1ull << 3, 0ull,
  0ull, 0ull, 1ull << 4, 0ull, 1ull << 5,
  0ull, 0ull, 0ull, 1ull << 6, 0ull,
  0ull, 0ull, 0ull, 0ull, 1ull << 7,
  0ull, 1ull << 8, 0ull, 0ull, 0ull,
  0ull, 0ull, 1ull << 9, 0ull, 0ull,
  0ull, 1ull << 10, 0ull, 1ull << 11, 0ull,
  0ull, 0ull, 1ull << 12, 0ull, 1ull << 13,
  0ull, 0ull, 0ull, 1ull << 14, 0ull,
  0ull, 0ull, 0ull, 0ull, 1ull << 15,
  0ull, 1ull << 16, 0ull, 0ull, 0ull,
  0ull, 0ull, 1ull << 17, 0ull, 0ull,
  0ull, 1ull << 18, 0ull, 1ull << 19, 0ull,
  0ull, 0ull, 1ull << 20, 0ull, 1ull << 21,
  0ull, 0ull, 0ull, 1ull << 22, 0ull,
  0ull, 0ull, 0ull, 0ull, 1ull << 23,
  0ull, 1ull << 24, 0ull, 0ull, 0ull,
  0ull, 0ull, 1ull << 25, 0ull, 0ull,
  0ull, 1ull << 26, 0ull, 1ull << 27, 0ull,
  0ull, 0ull, 1ull << 28, 0ull, 1ull << 29,
  0ull, 0ull, 0ull, 1ull << 30, 0ull,
  0ull, 0ull, 0ull, 0ull, 1ull << 31,
  0ull, 1ull << 32, 0ull, 0ull, 0ull,
  0ull, 0ull, 1ull << 33, 0ull, 0ull,
  0ull, 1ull << 34, 0ull, 1ull << 35, 0ull,
  0ull, 0ull, 1ull << 36, 0ull, 1ull << 37,
  0ull, 0ull, 0ull, 1ull << 38, 0ull,
  0ull, 0ull, 0ull, 0ull, 1ull << 39,
  0ull, 1ull << 40, 0ull, 0ull, 0ull,
  0ull, 0ull, 1ull << 41, 0ull, 0ull,
  0ull, 1ull << 42, 0ull, 1ull << 43, 0ull,
  0ull, 0ull, 1ull << 44, 0ull, 1ull << 45,
  0ull, 0ull, 0ull, 1ull << 46, 0ull,
  0ull, 0ull, 0ull, 0ull, 1ull << 47,
  0ull, 1ull << 48, 0ull, 0ull, 0ull,
  0ull, 0ull, 1ull << 49, 0ull, 0ull,
  0ull, 1ull << 50, 0ull, 1ull << 51, 0ull,
  0ull, 0ull, 1ull << 52, 0ull, 1ull << 53,
  0ull, 0ull, 0ull, 1ull << 54, 0ull,
  0ull, 0ull, 0ull, 0ull, 1ull << 55,
  0ull, 1ull << 56, 0ull, 0ull, 0ull,
  0ull, 0ull, 1ull << 57, 0ull, 0ull,
  0ull, 1ull << 58, 0ull, 1ull << 59, 0ull,
  0ull, 0ull, 1ull << 60, 0ull, 1ull << 61,
  0ull, 0ull, 0ull, 1ull << 62, 0ull,
  0ull, 0ull, 0ull, 0ull, 1ull << 63
};

/// Bitmasks needed to unset a specific bit in the sieve array
const Array<uint64_t, 240> BitSieve240::unset_bit_ =
{
  ~0ull, ~(1ull << 0), ~0ull, ~0ull, ~0ull,
  ~0ull, ~0ull, ~(1ull << 1), ~0ull, ~0ull,
  ~0ull, ~(1ull << 2), ~0ull, ~(1ull << 3), ~0ull,
  ~0ull, ~0ull, ~(1ull << 4), ~0ull, ~(1ull << 5),
  ~0ull, ~0ull, ~0ull, ~(1ull << 6), ~0ull,
  ~0ull, ~0ull, ~0ull, ~0ull, ~(1ull << 7),
  ~0ull, ~(1ull << 8), ~0ull, ~0ull, ~0ull,
  ~0ull, ~0ull, ~(1ull << 9), ~0ull, ~0ull,
  ~0ull, ~(1ull << 10), ~0ull, ~(1ull << 11), ~0ull,
  ~0ull, ~0ull, ~(1ull << 12), ~0ull, ~(1ull << 13),
  ~0ull, ~0ull, ~0ull, ~(1ull << 14), ~0ull,
  ~0ull, ~0ull, ~0ull, ~0ull, ~(1ull << 15),
  ~0ull, ~(1ull << 16), ~0ull, ~0ull, ~0ull,
  ~0ull, ~0ull, ~(1ull << 17), ~0ull, ~0ull,
  ~0ull, ~(1ull << 18), ~0ull, ~(1ull << 19), ~0ull,
  ~0ull, ~0ull, ~(1ull << 20), ~0ull, ~(1ull << 21),
  ~0ull, ~0ull, ~0ull, ~(1ull << 22), ~0ull,
  ~0ull, ~0ull, ~0ull, ~0ull, ~(1ull << 23),
  ~0ull, ~(1ull << 24), ~0ull, ~0ull, ~0ull,
  ~0ull, ~0ull, ~(1ull << 25), ~0ull, ~0ull,
  ~0ull, ~(1ull << 26), ~0ull, ~(1ull << 27), ~0ull,
  ~0ull, ~0ull, ~(1ull << 28), ~0ull, ~(1ull << 29),
  ~0ull, ~0ull, ~0ull, ~(1ull << 30), ~0ull,
  ~0ull, ~0ull, ~0ull, ~0ull, ~(1ull << 31),
  ~0ull, ~(1ull << 32), ~0ull, ~0ull, ~0ull,
  ~0ull, ~0ull, ~(1ull << 33), ~0ull, ~0ull,
  ~0ull, ~(1ull << 34), ~0ull, ~(1ull << 35), ~0ull,
  ~0ull, ~0ull, ~(1ull << 36), ~0ull, ~(1ull << 37),
  ~0ull, ~0ull, ~0ull, ~(1ull << 38), ~0ull,
  ~0ull, ~0ull, ~0ull, ~0ull, ~(1ull << 39),
  ~0ull, ~(1ull << 40), ~0ull, ~0ull, ~0ull,
  ~0ull, ~0ull, ~(1ull << 41), ~0ull, ~0ull,
  ~0ull, ~(1ull << 42), ~0ull, ~(1ull << 43), ~0ull,
  ~0ull, ~0ull, ~(1ull << 44), ~0ull, ~(1ull << 45),
  ~0ull, ~0ull, ~0ull, ~(1ull << 46), ~0ull,
  ~0ull, ~0ull, ~0ull, ~0ull, ~(1ull << 47),
  ~0ull, ~(1ull << 48), ~0ull, ~0ull, ~0ull,
  ~0ull, ~0ull, ~(1ull << 49), ~0ull, ~0ull,
  ~0ull, ~(1ull << 50), ~0ull, ~(1ull << 51), ~0ull,
  ~0ull, ~0ull, ~(1ull << 52), ~0ull, ~(1ull << 53),
  ~0ull, ~0ull, ~0ull, ~(1ull << 54), ~0ull,
  ~0ull, ~0ull, ~0ull, ~0ull, ~(1ull << 55),
  ~0ull, ~(1ull << 56), ~0ull, ~0ull, ~0ull,
  ~0ull, ~0ull, ~(1ull << 57), ~0ull, ~0ull,
  ~0ull, ~(1ull << 58), ~0ull, ~(1ull << 59), ~0ull,
  ~0ull, ~0ull, ~(1ull << 60), ~0ull, ~(1ull << 61),
  ~0ull, ~0ull, ~0ull, ~(1ull << 62), ~0ull,
  ~0ull, ~0ull, ~0ull, ~0ull, ~(1ull << 63)
};

/// Unset bits < start
const Array<uint64_t, 240> BitSieve240::unset_smaller_ =
{
  unset_s(0), unset_s(1), unset_s(2), unset_s(3), unset_s(4),
  unset_s(5), unset_s(6), unset_s(7), unset_s(8), unset_s(9),
  unset_s(10), unset_s(11), unset_s(12), unset_s(13), unset_s(14),
  unset_s(15), unset_s(16), unset_s(17), unset_s(18), unset_s(19),
  unset_s(20), unset_s(21), unset_s(22), unset_s(23), unset_s(24),
  unset_s(25), unset_s(26), unset_s(27), unset_s(28), unset_s(29),
  unset_s(30), unset_s(31), unset_s(32), unset_s(33), unset_s(34),
  unset_s(35), unset_s(36), unset_s(37), unset_s(38), unset_s(39),
  unset_s(40), unset_s(41), unset_s(42), unset_s(43), unset_s(44),
  unset_s(45), unset_s(46), unset_s(47), unset_s(48), unset_s(49),
  unset_s(50), unset_s(51), unset_s(52), unset_s(53), unset_s(54),
  unset_s(55), unset_s(56), unset_s(57), unset_s(58), unset_s(59),
  unset_s(60), unset_s(61), unset_s(62), unset_s(63), unset_s(64),
  unset_s(65), unset_s(66), unset_s(67), unset_s(68), unset_s(69),
  unset_s(70), unset_s(71), unset_s(72), unset_s(73), unset_s(74),
  unset_s(75), unset_s(76), unset_s(77), unset_s(78), unset_s(79),
  unset_s(80), unset_s(81), unset_s(82), unset_s(83), unset_s(84),
  unset_s(85), unset_s(86), unset_s(87), unset_s(88), unset_s(89),
  unset_s(90), unset_s(91), unset_s(92), unset_s(93), unset_s(94),
  unset_s(95), unset_s(96), unset_s(97), unset_s(98), unset_s(99),
  unset_s(100), unset_s(101), unset_s(102), unset_s(103), unset_s(104),
  unset_s(105), unset_s(106), unset_s(107), unset_s(108), unset_s(109),
  unset_s(110), unset_s(111), unset_s(112), unset_s(113), unset_s(114),
  unset_s(115), unset_s(116), unset_s(117), unset_s(118), unset_s(119),
  unset_s(120), unset_s(121), unset_s(122), unset_s(123), unset_s(124),
  unset_s(125), unset_s(126), unset_s(127), unset_s(128), unset_s(129),
  unset_s(130), unset_s(131), unset_s(132), unset_s(133), unset_s(134),
  unset_s(135), unset_s(136), unset_s(137), unset_s(138), unset_s(139),
  unset_s(140), unset_s(141), unset_s(142), unset_s(143), unset_s(144),
  unset_s(145), unset_s(146), unset_s(147), unset_s(148), unset_s(149),
  unset_s(150), unset_s(151), unset_s(152), unset_s(153), unset_s(154),
  unset_s(155), unset_s(156), unset_s(157), unset_s(158), unset_s(159),
  unset_s(160), unset_s(161), unset_s(162), unset_s(163), unset_s(164),
  unset_s(165), unset_s(166), unset_s(167), unset_s(168), unset_s(169),
  unset_s(170), unset_s(171), unset_s(172), unset_s(173), unset_s(174),
  unset_s(175), unset_s(176), unset_s(177), unset_s(178), unset_s(179),
  unset_s(180), unset_s(181), unset_s(182), unset_s(183), unset_s(184),
  unset_s(185), unset_s(186), unset_s(187), unset_s(188), unset_s(189),
  unset_s(190), unset_s(191), unset_s(192), unset_s(193), unset_s(194),
  unset_s(195), unset_s(196), unset_s(197), unset_s(198), unset_s(199),
  unset_s(200), unset_s(201), unset_s(202), unset_s(203), unset_s(204),
  unset_s(205), unset_s(206), unset_s(207), unset_s(208), unset_s(209),
  unset_s(210), unset_s(211), unset_s(212), unset_s(213), unset_s(214),
  unset_s(215), unset_s(216), unset_s(217), unset_s(218), unset_s(219),
  unset_s(220), unset_s(221), unset_s(222), unset_s(223), unset_s(224),
  unset_s(225), unset_s(226), unset_s(227), unset_s(228), unset_s(229),
  unset_s(230), unset_s(231), unset_s(232), unset_s(233), unset_s(234),
  unset_s(235), unset_s(236), unset_s(237), unset_s(238), unset_s(239)
};

/// Unset bits > stop
const Array<uint64_t, 240> BitSieve240::unset_larger_ =
{
  unset_l(0), unset_l(1), unset_l(2), unset_l(3), unset_l(4),
  unset_l(5), unset_l(6), unset_l(7), unset_l(8), unset_l(9),
  unset_l(10), unset_l(11), unset_l(12), unset_l(13), unset_l(14),
  unset_l(15), unset_l(16), unset_l(17), unset_l(18), unset_l(19),
  unset_l(20), unset_l(21), unset_l(22), unset_l(23), unset_l(24),
  unset_l(25), unset_l(26), unset_l(27), unset_l(28), unset_l(29),
  unset_l(30), unset_l(31), unset_l(32), unset_l(33), unset_l(34),
  unset_l(35), unset_l(36), unset_l(37), unset_l(38), unset_l(39),
  unset_l(40), unset_l(41), unset_l(42), unset_l(43), unset_l(44),
  unset_l(45), unset_l(46), unset_l(47), unset_l(48), unset_l(49),
  unset_l(50), unset_l(51), unset_l(52), unset_l(53), unset_l(54),
  unset_l(55), unset_l(56), unset_l(57), unset_l(58), unset_l(59),
  unset_l(60), unset_l(61), unset_l(62), unset_l(63), unset_l(64),
  unset_l(65), unset_l(66), unset_l(67), unset_l(68), unset_l(69),
  unset_l(70), unset_l(71), unset_l(72), unset_l(73), unset_l(74),
  unset_l(75), unset_l(76), unset_l(77), unset_l(78), unset_l(79),
  unset_l(80), unset_l(81), unset_l(82), unset_l(83), unset_l(84),
  unset_l(85), unset_l(86), unset_l(87), unset_l(88), unset_l(89),
  unset_l(90), unset_l(91), unset_l(92), unset_l(93), unset_l(94),
  unset_l(95), unset_l(96), unset_l(97), unset_l(98), unset_l(99),
  unset_l(100), unset_l(101), unset_l(102), unset_l(103), unset_l(104),
  unset_l(105), unset_l(106), unset_l(107), unset_l(108), unset_l(109),
  unset_l(110), unset_l(111), unset_l(112), unset_l(113), unset_l(114),
  unset_l(115), unset_l(116), unset_l(117), unset_l(118), unset_l(119),
  unset_l(120), unset_l(121), unset_l(122), unset_l(123), unset_l(124),
  unset_l(125), unset_l(126), unset_l(127), unset_l(128), unset_l(129),
  unset_l(130), unset_l(131), unset_l(132), unset_l(133), unset_l(134),
  unset_l(135), unset_l(136), unset_l(137), unset_l(138), unset_l(139),
  unset_l(140), unset_l(141), unset_l(142), unset_l(143), unset_l(144),
  unset_l(145), unset_l(146), unset_l(147), unset_l(148), unset_l(149),
  unset_l(150), unset_l(151), unset_l(152), unset_l(153), unset_l(154),
  unset_l(155), unset_l(156), unset_l(157), unset_l(158), unset_l(159),
  unset_l(160), unset_l(161), unset_l(162), unset_l(163), unset_l(164),
  unset_l(165), unset_l(166), unset_l(167), unset_l(168), unset_l(169),
  unset_l(170), unset_l(171), unset_l(172), unset_l(173), unset_l(174),
  unset_l(175), unset_l(176), unset_l(177), unset_l(178), unset_l(179),
  unset_l(180), unset_l(181), unset_l(182), unset_l(183), unset_l(184),
  unset_l(185), unset_l(186), unset_l(187), unset_l(188), unset_l(189),
  unset_l(190), unset_l(191), unset_l(192), unset_l(193), unset_l(194),
  unset_l(195), unset_l(196), unset_l(197), unset_l(198), unset_l(199),
  unset_l(200), unset_l(201), unset_l(202), unset_l(203), unset_l(204),
  unset_l(205), unset_l(206), unset_l(207), unset_l(208), unset_l(209),
  unset_l(210), unset_l(211), unset_l(212), unset_l(213), unset_l(214),
  unset_l(215), unset_l(216), unset_l(217), unset_l(218), unset_l(219),
  unset_l(220), unset_l(221), unset_l(222), unset_l(223), unset_l(224),
  unset_l(225), unset_l(226), unset_l(227), unset_l(228), unset_l(229),
  unset_l(230), unset_l(231), unset_l(232), unset_l(233), unset_l(234),
  unset_l(235), unset_l(236), unset_l(237), unset_l(238), unset_l(239)
};

} // namespace
