///
/// @file  Sieve.cpp
/// @brief This file implements a highly optimized prime sieving
///        algorithm for computing the special leaves (sometimes named
///        hard special leaves) in the combinatorial prime counting
///        algorithms (e.g. Lagarias-Miller-Odlyzko, Deleglise-Rivat,
///        Gourdon).
///
///        The Sieve class contains a sieve of Eratosthenes
///        implementation with 30 numbers per byte i.e. the 8 bits of
///        each byte correspond to the offsets: { 1, 7, 11, 13, 17,
///        19, 23, 29 }. Unlike a traditional prime sieve this sieve
///        is designed for use in the combinatorial prime counting
///        algorithms: this sieve removes primes as well as multiples
///        of primes and it counts the number of elements that have
///        been crossed off for the first time in the sieve array.
///
///        Since there is a large number of leaves for which we have
///        to count the number of unsieved elements in the sieve
///        array, Lagarias-Miller-Odlyzko have suggested using a
///        binary indexed tree data structure (a.k.a. Fenwick tree) to
///        speedup counting. However using a binary indexed tree is
///        bad for performance as it causes many cache misses and
///        branch mispredictions. For this reason this implementation
///        instead uses a linear counter array whose elements contain
///        the total count of unsieved elements in a certain interval.
///
///        In-depth description of this algorithm:
///        https://github.com/kimwalisch/primecount/blob/master/doc/Hard-Special-Leaves.md
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <Sieve.hpp>
#include <cpu_arch_macros.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <macros.hpp>
#include <min.hpp>
#include <Vector.hpp>
#include <popcnt.hpp>

#include <stdint.h>
#include <algorithm>
#include <cmath>

#if defined(ENABLE_ARM_SVE)
  #include <arm_sve.h>
#elif defined(ENABLE_AVX512_VPOPCNT)
  #include <immintrin.h>
#elif defined(ENABLE_MULTIARCH_ARM_SVE)
  #include <arm_sve.h>
  #include <cpu_supports_arm_sve.hpp>
#elif defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  #include <immintrin.h>
  #include <cpu_supports_avx512_vpopcnt.hpp>
#endif

namespace {

struct WheelInit
{
  uint8_t factor;
  uint8_t index;
};

using primecount::Array;

/// Categorize sieving primes according to their modulo 30
/// congruence class { 1, 7, 11, 13, 17, 19, 23, 29 }.
///
const Array<uint8_t, 30> wheel_offsets =
{
  0, 8 * 0, 0, 0, 0, 0,
  0, 8 * 1, 0, 0, 0, 8 * 2,
  0, 8 * 3, 0, 0, 0, 8 * 4,
  0, 8 * 5, 0, 0, 0, 8 * 6,
  0, 0,     0, 0, 0, 8 * 7
};

/// Used to calculate the first multiple > start of a
/// sieving prime that is coprime to 2, 3, 5.
///
const Array<WheelInit, 30> wheel_init
{{
  {1,  0}, {0,  0}, {5,  1}, {4,  1}, {3,  1},
  {2,  1}, {1,  1}, {0,  1}, {3,  2}, {2,  2},
  {1,  2}, {0,  2}, {1,  3}, {0,  3}, {3,  4},
  {2,  4}, {1,  4}, {0,  4}, {1,  5}, {0,  5},
  {3,  6}, {2,  6}, {1,  6}, {0,  6}, {5,  7},
  {4,  7}, {3,  7}, {2,  7}, {1,  7}, {0,  7}
}};

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

#if defined(__BYTE_ORDER__) && \
    defined(__ORDER_BIG_ENDIAN__) && \
    __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__

/// Since the sieve array is a byte array but the
/// unset_smaller and unset_larger arrays are of
/// type uint64_t we need to reverse the byte order
/// of the unset_smaller and unset_larger arrays
/// on big endian CPU architectures.

/// Reverse byte order of 64-bit integer
constexpr uint64_t bswap64(uint64_t i, uint64_t j = 0, uint64_t n = 0)
{
  return (n == sizeof(uint64_t)) ? j :
    bswap64(i >> 8, (j << 8) | (i & 0xff), n + 1);
}

/// Bitmask to unset bits <= n
constexpr uint64_t unset_s(int n)
{
  return bswap64(~0ull << left_shift(n));
}

/// Bitmask to unset bits >= n
constexpr uint64_t unset_l(int n)
{
  return bswap64((n == 0) ? 0 : ~0ull >> right_shift(n));
}

#else

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

#endif

} // namespace

namespace primecount {

/// Unset bits < start
const Array<uint64_t, 240> Sieve::unset_smaller =
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
const Array<uint64_t, 240> Sieve::unset_larger =
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

Sieve::Sieve(uint64_t low,
             uint64_t segment_size, 
             uint64_t wheel_size)
{
  ASSERT(low % 30 == 0);
  ASSERT(segment_size % 240 == 0);

  start_ = low;
  segment_size = align_segment_size(segment_size);

  // sieve_size = segment_size / 30 as each byte corresponds
  // to 30 numbers i.e. the 8 bits correspond to the
  // offsets = {1, 7, 11, 13, 17, 19, 23, 29}.
  sieve_.resize(segment_size / 30);
  wheel_.reserve(wheel_size);
  wheel_.resize(4);
  allocate_counter(low);
}

/// Each element of the counter array contains the current
/// number of unsieved elements in the interval:
/// [i * counter_.dist, (i + 1) * counter_.dist[.
/// Ideally each element of the counter array should
/// represent an interval of size O(sqrt(average_leaf_dist)).
/// Also the counter distance should be adjusted regularly
/// whilst sieving as the distance between consecutive
/// leaves is very small ~ log(x) at the beginning of the
/// sieving algorithm but grows up to segment_size towards
/// the end of the algorithm.
///
void Sieve::allocate_counter(uint64_t low)
{
  // Default count_popcnt64() algorithm
  uint64_t bytes_per_instruction = sizeof(uint64_t);

#if defined(ENABLE_AVX512_VPOPCNT)
  // count_avx512() algorithm
  bytes_per_instruction = sizeof(__m512i);
#elif defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  // count_avx512() algorithm
  if (cpu_supports_avx512_vpopcnt)
    bytes_per_instruction = sizeof(__m512i);
#elif defined(ENABLE_ARM_SVE)
  // count_arm_sve() algorithm
  bytes_per_instruction = svcntd() * sizeof(uint64_t);
#elif defined(ENABLE_MULTIARCH_ARM_SVE)
  // count_arm_sve() algorithm
  if (cpu_supports_sve)
    bytes_per_instruction = svcntd() * sizeof(uint64_t);
#endif

  double average_leaf_dist = std::sqrt(low);
  double counter_dist = std::sqrt(average_leaf_dist);

  // Here we balance counting with the counter array and
  // counting from the sieve array using the 64-bit POPCNT
  // instruction. Since the 64-bit POPCNT instructions
  // allows to count a distance of 240 using a single
  // instruction we slightly increase the counter distance
  // and slightly decrease the size of the counter array.
  ASSERT(bytes_per_instruction >= sizeof(uint64_t));
  uint64_t dist_per_instruction = bytes_per_instruction * 30;
  counter_.dist = uint64_t(counter_dist * std::sqrt(dist_per_instruction));

  // Increasing the minimum counter distance decreases the
  // branch mispredictions (good) but on the other hand
  // increases the number of executed instructions (bad).
  // In my benchmarks setting the minimum amount of bytes to
  // bytes_per_instruction * 8 (or 16) performed best.
  uint64_t bytes = counter_.dist / 30;
  bytes = max(bytes, bytes_per_instruction * 8);
  bytes = next_power_of_2(bytes);

  // Make sure the counter (32-bit) doesn't overflow.
  // This can never happen since each counter array element
  // only counts the number of unsieved elements (1 bits) in
  // an interval of size: sieve_limit^(1/4) * sqrt(240).
  // Hence the max(counter value) = 2^18.
  ASSERT(bytes * 8 <= pstd::numeric_limits<uint32_t>::max());
  uint64_t counter_size = ceil_div(sieve_.size(), bytes);
  counter_.counter.resize(counter_size);
  counter_.dist = bytes * 30;
  counter_.log2_dist = ilog2(bytes);
}

/// The segment size is sieve.size() * 30 as each
/// byte corresponds to 30 numbers.
///
uint64_t Sieve::segment_size() const
{
  return sieve_.size() * 30;
}

/// segment_size must be a multiple of 240 as we
/// process 64-bit words (8 bytes) and each
/// byte contains 30 numbers.
///
uint64_t Sieve::align_segment_size(uint64_t size)
{
  size = max(size, 240);

  if (size % 240)
    size += 240 - size % 240;

  return size;
}

void Sieve::reset_sieve(uint64_t low, uint64_t high)
{
  std::fill_n(sieve_.data(), sieve_.size(), 0xff);
  uint64_t size = high - low;

  if (size < segment_size())
  {
    uint64_t last = size - 1;
    size = align_segment_size(size);
    sieve_.resize(size / 30);
    auto sieve64 = (uint64_t*) sieve_.data();
    sieve64[last / 240] &= unset_larger[last % 240];
  }
}

void Sieve::reset_counter()
{
  prev_stop_ = 0;
  count_ = 0;
  counter_.i = 0;
  counter_.sum = 0;
  counter_.stop = counter_.dist;
}

void Sieve::init_counter(uint64_t low, uint64_t high)
{
  reset_counter();
  total_count_ = 0;

  uint64_t start = 0;
  uint64_t max_stop = (high - 1) - low;

  while (start <= max_stop)
  {
    uint64_t stop = start + counter_.dist - 1;
    stop = min(stop, max_stop);
    uint64_t cnt = count(start, stop);
    uint64_t byte_index = start / 30;
    uint64_t i = byte_index >> counter_.log2_dist;

    counter_[i] = (uint32_t) cnt;
    total_count_ += cnt;
    start += counter_.dist;
  }
}

/// Count 1 bits inside [start, stop]
uint64_t Sieve::count(uint64_t start, uint64_t stop) const
{
  #if defined(ENABLE_ARM_SVE)
    return count_arm_sve(start, stop);
  #elif defined(ENABLE_AVX512_VPOPCNT)
    return count_avx512(start, stop);
  #elif defined(ENABLE_MULTIARCH_ARM_SVE)
    return cpu_supports_sve ? count_arm_sve(start, stop) : count_popcnt64(start, stop);
  #elif defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
    return cpu_supports_avx512_vpopcnt ? count_avx512(start, stop) : count_popcnt64(start, stop);
  #else
    return count_popcnt64(start, stop);
  #endif
}

#if defined(ENABLE_PORTABLE_POPCNT64)

/// Count 1 bits inside [start, stop].
/// The distance [start, stop] is small here < sqrt(segment_size),
/// hence we simply count the number of unsieved elements
/// by linearly iterating over the sieve array.
///
uint64_t Sieve::count_popcnt64(uint64_t start, uint64_t stop) const
{
  if (start > stop)
    return 0;

  ASSERT(stop - start < segment_size());
  uint64_t start_idx = start / 240;
  uint64_t stop_idx = stop / 240;
  uint64_t m1 = unset_smaller[start % 240];
  uint64_t m2 = unset_larger[stop % 240];

  // Branchfree bitmask calculation:
  // m1 = (start_idx != stop_idx) ? m1 : m1 & m2;
  m1 &= (-(start_idx != stop_idx) | m2);
  // m2 = (start_idx != stop_idx) ? m2 : 0;
  m2 &= -(start_idx != stop_idx);

  const uint64_t* sieve64 = (const uint64_t*) sieve_.data();
  uint64_t start_bits = sieve64[start_idx] & m1;
  uint64_t stop_bits = sieve64[stop_idx] & m2;
  uint64_t cnt = popcnt64(start_bits);
  cnt += popcnt64(stop_bits);

  for (uint64_t i = start_idx + 1; i < stop_idx; i++)
    cnt += popcnt64(sieve64[i]);

  return cnt;
}

#endif

#if defined(ENABLE_AVX512_VPOPCNT) || \
    defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)

/// Count 1 bits inside [start, stop].
/// The distance [start, stop] is small here < sqrt(segment_size),
/// hence we simply count the number of unsieved elements
/// by linearly iterating over the sieve array.
///
#if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  __attribute__ ((target ("avx512f,avx512vpopcntdq")))
#endif
uint64_t Sieve::count_avx512(uint64_t start, uint64_t stop) const
{
  if (start > stop)
    return 0;

  ASSERT(stop - start < segment_size());
  uint64_t start_idx = start / 240;
  uint64_t stop_idx = stop / 240;
  uint64_t m1 = unset_smaller[start % 240];
  uint64_t m2 = unset_larger[stop % 240];

  // Branchfree bitmask calculation:
  // m1 = (start_idx != stop_idx) ? m1 : m1 & m2;
  m1 &= (-(start_idx != stop_idx) | m2);
  // m2 = (start_idx != stop_idx) ? m2 : 0;
  m2 &= -(start_idx != stop_idx);

  const uint64_t* sieve64 = (const uint64_t*) sieve_.data();
  uint64_t start_bits = sieve64[start_idx] & m1;
  uint64_t stop_bits = sieve64[stop_idx] & m2;
  __m512i vec = _mm512_set_epi64(0, 0, 0, 0, 0, 0, stop_bits, start_bits);
  __m512i vcnt = _mm512_popcnt_epi64(vec);
  uint64_t i = start_idx + 1;

  // Compute this for loop using AVX512.
  // for (i = start_idx + 1; i < stop_idx; i++)
  //   cnt += popcnt64(sieve64[i]);

  for (; i + 8 < stop_idx; i += 8)
  {
    vec = _mm512_loadu_epi64(&sieve64[i]);
    vec = _mm512_popcnt_epi64(vec);
    vcnt = _mm512_add_epi64(vcnt, vec);
  }

  __mmask8 mask = (__mmask8) (0xff >> (i + 8 - stop_idx));
  vec = _mm512_maskz_loadu_epi64(mask, &sieve64[i]);
  vec = _mm512_popcnt_epi64(vec);
  vcnt = _mm512_add_epi64(vcnt, vec);
  return _mm512_reduce_add_epi64(vcnt);
}

#elif defined(ENABLE_ARM_SVE) || \
      defined(ENABLE_MULTIARCH_ARM_SVE)

/// Count 1 bits inside [start, stop].
/// The distance [start, stop] is small here < sqrt(segment_size),
/// hence we simply count the number of unsieved elements
/// by linearly iterating over the sieve array.
///
#if defined(ENABLE_MULTIARCH_ARM_SVE)
  __attribute__ ((target ("arch=armv8-a+sve")))
#endif
uint64_t Sieve::count_arm_sve(uint64_t start, uint64_t stop) const
{
  if (start > stop)
    return 0;

  ASSERT(stop - start < segment_size());
  uint64_t start_idx = start / 240;
  uint64_t stop_idx = stop / 240;
  uint64_t m1 = unset_smaller[start % 240];
  uint64_t m2 = unset_larger[stop % 240];

  // Branchfree bitmask calculation:
  // m1 = (start_idx != stop_idx) ? m1 : m1 & m2;
  m1 &= (-(start_idx != stop_idx) | m2);
  // m2 = (start_idx != stop_idx) ? m2 : 0;
  m2 &= -(start_idx != stop_idx);

  const uint64_t* sieve64 = (const uint64_t*) sieve_.data();
  uint64_t start_bits = sieve64[start_idx] & m1;
  uint64_t stop_bits = sieve64[stop_idx] & m2;
  ASSERT(svcntd() >= 2);
  svuint64_t vec = svinsr_n_u64(svdup_u64(start_bits), stop_bits);
  svuint64_t vcnt = svcnt_u64_z(svwhilelt_b64(0, 2), vec);
  uint64_t i = start_idx + 1;

  // Compute this for loop using ARM SVE.
  // for (i = start_idx + 1; i < stop_idx; i++)
  //   cnt += popcnt64(sieve64[i]);

  for (; i + svcntd() < stop_idx; i += svcntd())
  {
    vec = svld1_u64(svptrue_b64(), &sieve64[i]);
    vec = svcnt_u64_x(svptrue_b64(), vec);
    vcnt = svadd_u64_x(svptrue_b64(), vcnt, vec);
  }

  svbool_t pg = svwhilelt_b64(i, stop_idx);
  vec = svld1_u64(pg, &sieve64[i]);
  vec = svcnt_u64_z(pg, vec);
  vcnt = svadd_u64_x(svptrue_b64(), vcnt, vec);
  return svaddv_u64(svptrue_b64(), vcnt);
}

#endif

/// Add a sieving prime to the sieve.
/// Calculates the first multiple > start of prime that
/// is not divisible by 2, 3, 5 and its wheel index.
///
void Sieve::add(uint64_t prime)
{
  ASSERT(start_ % 30 == 0);

  // first multiple > start_
  uint64_t quotient = start_ / prime + 1;
  uint64_t multiple = prime * quotient;

  // find next multiple of prime that
  // is not divisible by 2, 3, 5
  uint64_t factor = wheel_init[quotient % 30].factor;
  multiple += prime * factor;

  ASSERT(multiple % 2 != 0);
  ASSERT(multiple % 3 != 0);
  ASSERT(multiple % 5 != 0);

  multiple = (multiple - start_) / 30;
  uint32_t multiple32 = (uint32_t) multiple;

  // calculate wheel index of multiple
  uint32_t index = wheel_init[quotient % 30].index;
  index += wheel_offsets[prime % 30];
  wheel_.emplace_back(multiple32, index);
}

/// Remove the i-th prime and the multiples of the i-th prime
/// from the sieve array. Used for pre-sieving.
///
void Sieve::cross_off(uint64_t prime, uint64_t i)
{
  if (i >= wheel_.size())
    add(prime);

  prime /= 30;
  Wheel& wheel = wheel_[i];
  uint64_t m = wheel.multiple;
  uint8_t* sieve = sieve_.data();
  uint64_t sieve_size = sieve_.size();

  #define CHECK_FINISHED(wheel_index) \
    if_unlikely(m >= sieve_size) \
    { \
      wheel.index = wheel_index; \
      wheel.multiple = (uint32_t) (m - sieve_size); \
      return; \
    }

  ASSERT(wheel.index <= 63);

  switch (wheel.index)
  {
    for (;;)
    {
      case 0: {
                uint64_t max_offset = m + prime * 28;
                uint64_t limit = std::max(max_offset, sieve_size) - max_offset;

                for (; m < limit; m += prime * 30 + 1)
                {
                  sieve[m + prime *  0] &= ~(1 << 0);
                  sieve[m + prime *  6] &= ~(1 << 1);
                  sieve[m + prime * 10] &= ~(1 << 2);
                  sieve[m + prime * 12] &= ~(1 << 3);
                  sieve[m + prime * 16] &= ~(1 << 4);
                  sieve[m + prime * 18] &= ~(1 << 5);
                  sieve[m + prime * 22] &= ~(1 << 6);
                  sieve[m + prime * 28] &= ~(1 << 7);
                }
              }
              CHECK_FINISHED(0); sieve[m] &= ~(1 << 0); m += prime * 6 + 0; FALLTHROUGH;
      case 1: CHECK_FINISHED(1); sieve[m] &= ~(1 << 1); m += prime * 4 + 0; FALLTHROUGH;
      case 2: CHECK_FINISHED(2); sieve[m] &= ~(1 << 2); m += prime * 2 + 0; FALLTHROUGH;
      case 3: CHECK_FINISHED(3); sieve[m] &= ~(1 << 3); m += prime * 4 + 0; FALLTHROUGH;
      case 4: CHECK_FINISHED(4); sieve[m] &= ~(1 << 4); m += prime * 2 + 0; FALLTHROUGH;
      case 5: CHECK_FINISHED(5); sieve[m] &= ~(1 << 5); m += prime * 4 + 0; FALLTHROUGH;
      case 6: CHECK_FINISHED(6); sieve[m] &= ~(1 << 6); m += prime * 6 + 0; FALLTHROUGH;
      case 7: CHECK_FINISHED(7); sieve[m] &= ~(1 << 7); m += prime * 2 + 1;
    }

    for (;;)
    {
      case  8: {
                 uint64_t max_offset = m + prime * 28 + 6;
                 uint64_t limit = std::max(max_offset, sieve_size) - max_offset;

                 for (; m < limit; m += prime * 30 + 7)
                 {
                   sieve[m + prime *  0 + 0] &= ~(1 << 1);
                   sieve[m + prime *  6 + 1] &= ~(1 << 5);
                   sieve[m + prime * 10 + 2] &= ~(1 << 4);
                   sieve[m + prime * 12 + 3] &= ~(1 << 0);
                   sieve[m + prime * 16 + 3] &= ~(1 << 7);
                   sieve[m + prime * 18 + 4] &= ~(1 << 3);
                   sieve[m + prime * 22 + 5] &= ~(1 << 2);
                   sieve[m + prime * 28 + 6] &= ~(1 << 6);
                 }
               }
               CHECK_FINISHED( 8); sieve[m] &= ~(1 << 1); m += prime * 6 + 1; FALLTHROUGH;
      case  9: CHECK_FINISHED( 9); sieve[m] &= ~(1 << 5); m += prime * 4 + 1; FALLTHROUGH;
      case 10: CHECK_FINISHED(10); sieve[m] &= ~(1 << 4); m += prime * 2 + 1; FALLTHROUGH;
      case 11: CHECK_FINISHED(11); sieve[m] &= ~(1 << 0); m += prime * 4 + 0; FALLTHROUGH;
      case 12: CHECK_FINISHED(12); sieve[m] &= ~(1 << 7); m += prime * 2 + 1; FALLTHROUGH;
      case 13: CHECK_FINISHED(13); sieve[m] &= ~(1 << 3); m += prime * 4 + 1; FALLTHROUGH;
      case 14: CHECK_FINISHED(14); sieve[m] &= ~(1 << 2); m += prime * 6 + 1; FALLTHROUGH;
      case 15: CHECK_FINISHED(15); sieve[m] &= ~(1 << 6); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 16: {
                 uint64_t max_offset = m + prime * 28 + 10;
                 uint64_t limit = std::max(max_offset, sieve_size) - max_offset;

                 for (; m < limit; m += prime * 30 + 11)
                 {
                   sieve[m + prime *  0 +  0] &= ~(1 << 2);
                   sieve[m + prime *  6 +  2] &= ~(1 << 4);
                   sieve[m + prime * 10 +  4] &= ~(1 << 0);
                   sieve[m + prime * 12 +  4] &= ~(1 << 6);
                   sieve[m + prime * 16 +  6] &= ~(1 << 1);
                   sieve[m + prime * 18 +  6] &= ~(1 << 7);
                   sieve[m + prime * 22 +  8] &= ~(1 << 3);
                   sieve[m + prime * 28 + 10] &= ~(1 << 5);
                 }
               }
               CHECK_FINISHED(16); sieve[m] &= ~(1 << 2); m += prime * 6 + 2; FALLTHROUGH;
      case 17: CHECK_FINISHED(17); sieve[m] &= ~(1 << 4); m += prime * 4 + 2; FALLTHROUGH;
      case 18: CHECK_FINISHED(18); sieve[m] &= ~(1 << 0); m += prime * 2 + 0; FALLTHROUGH;
      case 19: CHECK_FINISHED(19); sieve[m] &= ~(1 << 6); m += prime * 4 + 2; FALLTHROUGH;
      case 20: CHECK_FINISHED(20); sieve[m] &= ~(1 << 1); m += prime * 2 + 0; FALLTHROUGH;
      case 21: CHECK_FINISHED(21); sieve[m] &= ~(1 << 7); m += prime * 4 + 2; FALLTHROUGH;
      case 22: CHECK_FINISHED(22); sieve[m] &= ~(1 << 3); m += prime * 6 + 2; FALLTHROUGH;
      case 23: CHECK_FINISHED(23); sieve[m] &= ~(1 << 5); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 24: {
                 uint64_t max_offset = m + prime * 28 + 12;
                 uint64_t limit = std::max(max_offset, sieve_size) - max_offset;

                 for (; m < limit; m += prime * 30 + 13)
                 {
                   sieve[m + prime *  0 +  0] &= ~(1 << 3);
                   sieve[m + prime *  6 +  3] &= ~(1 << 0);
                   sieve[m + prime * 10 +  4] &= ~(1 << 6);
                   sieve[m + prime * 12 +  5] &= ~(1 << 5);
                   sieve[m + prime * 16 +  7] &= ~(1 << 2);
                   sieve[m + prime * 18 +  8] &= ~(1 << 1);
                   sieve[m + prime * 22 +  9] &= ~(1 << 7);
                   sieve[m + prime * 28 + 12] &= ~(1 << 4);
                 }
               }
               CHECK_FINISHED(24); sieve[m] &= ~(1 << 3); m += prime * 6 + 3; FALLTHROUGH;
      case 25: CHECK_FINISHED(25); sieve[m] &= ~(1 << 0); m += prime * 4 + 1; FALLTHROUGH;
      case 26: CHECK_FINISHED(26); sieve[m] &= ~(1 << 6); m += prime * 2 + 1; FALLTHROUGH;
      case 27: CHECK_FINISHED(27); sieve[m] &= ~(1 << 5); m += prime * 4 + 2; FALLTHROUGH;
      case 28: CHECK_FINISHED(28); sieve[m] &= ~(1 << 2); m += prime * 2 + 1; FALLTHROUGH;
      case 29: CHECK_FINISHED(29); sieve[m] &= ~(1 << 1); m += prime * 4 + 1; FALLTHROUGH;
      case 30: CHECK_FINISHED(30); sieve[m] &= ~(1 << 7); m += prime * 6 + 3; FALLTHROUGH;
      case 31: CHECK_FINISHED(31); sieve[m] &= ~(1 << 4); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 32: {
                 uint64_t max_offset = m + prime * 28 + 16;
                 uint64_t limit = std::max(max_offset, sieve_size) - max_offset;

                 for (; m < limit; m += prime * 30 + 17)
                 {
                   sieve[m + prime *  0 +  0] &= ~(1 << 4);
                   sieve[m + prime *  6 +  3] &= ~(1 << 7);
                   sieve[m + prime * 10 +  6] &= ~(1 << 1);
                   sieve[m + prime * 12 +  7] &= ~(1 << 2);
                   sieve[m + prime * 16 +  9] &= ~(1 << 5);
                   sieve[m + prime * 18 + 10] &= ~(1 << 6);
                   sieve[m + prime * 22 + 13] &= ~(1 << 0);
                   sieve[m + prime * 28 + 16] &= ~(1 << 3);
                 }
               }
               CHECK_FINISHED(32); sieve[m] &= ~(1 << 4); m += prime * 6 + 3; FALLTHROUGH;
      case 33: CHECK_FINISHED(33); sieve[m] &= ~(1 << 7); m += prime * 4 + 3; FALLTHROUGH;
      case 34: CHECK_FINISHED(34); sieve[m] &= ~(1 << 1); m += prime * 2 + 1; FALLTHROUGH;
      case 35: CHECK_FINISHED(35); sieve[m] &= ~(1 << 2); m += prime * 4 + 2; FALLTHROUGH;
      case 36: CHECK_FINISHED(36); sieve[m] &= ~(1 << 5); m += prime * 2 + 1; FALLTHROUGH;
      case 37: CHECK_FINISHED(37); sieve[m] &= ~(1 << 6); m += prime * 4 + 3; FALLTHROUGH;
      case 38: CHECK_FINISHED(38); sieve[m] &= ~(1 << 0); m += prime * 6 + 3; FALLTHROUGH;
      case 39: CHECK_FINISHED(39); sieve[m] &= ~(1 << 3); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 40: {
                 uint64_t max_offset = m + prime * 28 + 18;
                 uint64_t limit = std::max(max_offset, sieve_size) - max_offset;

                 for (; m < limit; m += prime * 30 + 19)
                 {
                   sieve[m + prime *  0 +  0] &= ~(1 << 5);
                   sieve[m + prime *  6 +  4] &= ~(1 << 3);
                   sieve[m + prime * 10 +  6] &= ~(1 << 7);
                   sieve[m + prime * 12 +  8] &= ~(1 << 1);
                   sieve[m + prime * 16 + 10] &= ~(1 << 6);
                   sieve[m + prime * 18 + 12] &= ~(1 << 0);
                   sieve[m + prime * 22 + 14] &= ~(1 << 4);
                   sieve[m + prime * 28 + 18] &= ~(1 << 2);
                 }
               }
               CHECK_FINISHED(40); sieve[m] &= ~(1 << 5); m += prime * 6 + 4; FALLTHROUGH;
      case 41: CHECK_FINISHED(41); sieve[m] &= ~(1 << 3); m += prime * 4 + 2; FALLTHROUGH;
      case 42: CHECK_FINISHED(42); sieve[m] &= ~(1 << 7); m += prime * 2 + 2; FALLTHROUGH;
      case 43: CHECK_FINISHED(43); sieve[m] &= ~(1 << 1); m += prime * 4 + 2; FALLTHROUGH;
      case 44: CHECK_FINISHED(44); sieve[m] &= ~(1 << 6); m += prime * 2 + 2; FALLTHROUGH;
      case 45: CHECK_FINISHED(45); sieve[m] &= ~(1 << 0); m += prime * 4 + 2; FALLTHROUGH;
      case 46: CHECK_FINISHED(46); sieve[m] &= ~(1 << 4); m += prime * 6 + 4; FALLTHROUGH;
      case 47: CHECK_FINISHED(47); sieve[m] &= ~(1 << 2); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 48: {
                 uint64_t max_offset = m + prime * 28 + 22;
                 uint64_t limit = std::max(max_offset, sieve_size) - max_offset;

                 for (; m < limit; m += prime * 30 + 23)
                 {
                   sieve[m + prime *  0 +  0] &= ~(1 << 6);
                   sieve[m + prime *  6 +  5] &= ~(1 << 2);
                   sieve[m + prime * 10 +  8] &= ~(1 << 3);
                   sieve[m + prime * 12 +  9] &= ~(1 << 7);
                   sieve[m + prime * 16 + 13] &= ~(1 << 0);
                   sieve[m + prime * 18 + 14] &= ~(1 << 4);
                   sieve[m + prime * 22 + 17] &= ~(1 << 5);
                   sieve[m + prime * 28 + 22] &= ~(1 << 1);
                 }
               }
               CHECK_FINISHED(48); sieve[m] &= ~(1 << 6); m += prime * 6 + 5; FALLTHROUGH;
      case 49: CHECK_FINISHED(49); sieve[m] &= ~(1 << 2); m += prime * 4 + 3; FALLTHROUGH;
      case 50: CHECK_FINISHED(50); sieve[m] &= ~(1 << 3); m += prime * 2 + 1; FALLTHROUGH;
      case 51: CHECK_FINISHED(51); sieve[m] &= ~(1 << 7); m += prime * 4 + 4; FALLTHROUGH;
      case 52: CHECK_FINISHED(52); sieve[m] &= ~(1 << 0); m += prime * 2 + 1; FALLTHROUGH;
      case 53: CHECK_FINISHED(53); sieve[m] &= ~(1 << 4); m += prime * 4 + 3; FALLTHROUGH;
      case 54: CHECK_FINISHED(54); sieve[m] &= ~(1 << 5); m += prime * 6 + 5; FALLTHROUGH;
      case 55: CHECK_FINISHED(55); sieve[m] &= ~(1 << 1); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 56: {
                 uint64_t max_offset = m + prime * 28 + 28;
                 uint64_t limit = std::max(max_offset, sieve_size) - max_offset;

                 for (; m < limit; m += prime * 30 + 29)
                 {
                   sieve[m + prime *  0 +  0] &= ~(1 << 7);
                   sieve[m + prime *  6 +  6] &= ~(1 << 6);
                   sieve[m + prime * 10 + 10] &= ~(1 << 5);
                   sieve[m + prime * 12 + 12] &= ~(1 << 4);
                   sieve[m + prime * 16 + 16] &= ~(1 << 3);
                   sieve[m + prime * 18 + 18] &= ~(1 << 2);
                   sieve[m + prime * 22 + 22] &= ~(1 << 1);
                   sieve[m + prime * 28 + 28] &= ~(1 << 0);
                 }
               }
               CHECK_FINISHED(56); sieve[m] &= ~(1 << 7); m += prime * 6 + 6; FALLTHROUGH;
      case 57: CHECK_FINISHED(57); sieve[m] &= ~(1 << 6); m += prime * 4 + 4; FALLTHROUGH;
      case 58: CHECK_FINISHED(58); sieve[m] &= ~(1 << 5); m += prime * 2 + 2; FALLTHROUGH;
      case 59: CHECK_FINISHED(59); sieve[m] &= ~(1 << 4); m += prime * 4 + 4; FALLTHROUGH;
      case 60: CHECK_FINISHED(60); sieve[m] &= ~(1 << 3); m += prime * 2 + 2; FALLTHROUGH;
      case 61: CHECK_FINISHED(61); sieve[m] &= ~(1 << 2); m += prime * 4 + 4; FALLTHROUGH;
      case 62: CHECK_FINISHED(62); sieve[m] &= ~(1 << 1); m += prime * 6 + 6; FALLTHROUGH;
      case 63: CHECK_FINISHED(63); sieve[m] &= ~(1 << 0); m += prime * 2 + 1;
    }

    default: UNREACHABLE;
  }

  #undef CHECK_FINISHED
}

/// Remove the i-th prime and the multiples of the i-th prime
/// from the sieve array. Also counts the number of elements
/// removed for the first time i.e. the count of sieved elements
/// whose least prime factor is the i-th prime.
///
void Sieve::cross_off_count(uint64_t prime, uint64_t i)
{
  if (i >= wheel_.size())
    add(prime);

  reset_counter();
  Wheel& wheel = wheel_[i];
  prime /= 30;

  uint64_t m = wheel.multiple;
  uint64_t total_count = total_count_;
  uint64_t counter_log2_dist = counter_.log2_dist;
  uint64_t sieve_size = sieve_.size();
  uint32_t* counter = &counter_[0];
  uint8_t* sieve = &sieve_[0];

  #define CHECK_FINISHED(wheel_index) \
    if_unlikely(m >= sieve_size) \
    { \
      wheel.index = wheel_index; \
      wheel.multiple = (uint32_t) (m - sieve_size); \
      total_count_ = total_count; \
      return; \
    }

  #define COUNT_UNSET_BIT(bit_index) \
    { \
      std::size_t sieve_byte = sieve[m]; \
      std::size_t is_bit = (sieve_byte >> bit_index) & 1; \
      sieve[m] &= ~(1 << bit_index); \
      counter[m >> counter_log2_dist] -= (uint32_t) is_bit; \
      total_count -= (uint64_t) is_bit; \
    }

  ASSERT(wheel.index <= 63);

  switch (wheel.index)
  {
    for (;;)
    {
      case 0: CHECK_FINISHED(0); COUNT_UNSET_BIT(0); m += prime * 6 + 0; FALLTHROUGH;
      case 1: CHECK_FINISHED(1); COUNT_UNSET_BIT(1); m += prime * 4 + 0; FALLTHROUGH;
      case 2: CHECK_FINISHED(2); COUNT_UNSET_BIT(2); m += prime * 2 + 0; FALLTHROUGH;
      case 3: CHECK_FINISHED(3); COUNT_UNSET_BIT(3); m += prime * 4 + 0; FALLTHROUGH;
      case 4: CHECK_FINISHED(4); COUNT_UNSET_BIT(4); m += prime * 2 + 0; FALLTHROUGH;
      case 5: CHECK_FINISHED(5); COUNT_UNSET_BIT(5); m += prime * 4 + 0; FALLTHROUGH;
      case 6: CHECK_FINISHED(6); COUNT_UNSET_BIT(6); m += prime * 6 + 0; FALLTHROUGH;
      case 7: CHECK_FINISHED(7); COUNT_UNSET_BIT(7); m += prime * 2 + 1;
    }

    for (;;)
    {
      case  8: CHECK_FINISHED( 8); COUNT_UNSET_BIT(1); m += prime * 6 + 1; FALLTHROUGH;
      case  9: CHECK_FINISHED( 9); COUNT_UNSET_BIT(5); m += prime * 4 + 1; FALLTHROUGH;
      case 10: CHECK_FINISHED(10); COUNT_UNSET_BIT(4); m += prime * 2 + 1; FALLTHROUGH;
      case 11: CHECK_FINISHED(11); COUNT_UNSET_BIT(0); m += prime * 4 + 0; FALLTHROUGH;
      case 12: CHECK_FINISHED(12); COUNT_UNSET_BIT(7); m += prime * 2 + 1; FALLTHROUGH;
      case 13: CHECK_FINISHED(13); COUNT_UNSET_BIT(3); m += prime * 4 + 1; FALLTHROUGH;
      case 14: CHECK_FINISHED(14); COUNT_UNSET_BIT(2); m += prime * 6 + 1; FALLTHROUGH;
      case 15: CHECK_FINISHED(15); COUNT_UNSET_BIT(6); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 16: CHECK_FINISHED(16); COUNT_UNSET_BIT(2); m += prime * 6 + 2; FALLTHROUGH;
      case 17: CHECK_FINISHED(17); COUNT_UNSET_BIT(4); m += prime * 4 + 2; FALLTHROUGH;
      case 18: CHECK_FINISHED(18); COUNT_UNSET_BIT(0); m += prime * 2 + 0; FALLTHROUGH;
      case 19: CHECK_FINISHED(19); COUNT_UNSET_BIT(6); m += prime * 4 + 2; FALLTHROUGH;
      case 20: CHECK_FINISHED(20); COUNT_UNSET_BIT(1); m += prime * 2 + 0; FALLTHROUGH;
      case 21: CHECK_FINISHED(21); COUNT_UNSET_BIT(7); m += prime * 4 + 2; FALLTHROUGH;
      case 22: CHECK_FINISHED(22); COUNT_UNSET_BIT(3); m += prime * 6 + 2; FALLTHROUGH;
      case 23: CHECK_FINISHED(23); COUNT_UNSET_BIT(5); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 24: CHECK_FINISHED(24); COUNT_UNSET_BIT(3); m += prime * 6 + 3; FALLTHROUGH;
      case 25: CHECK_FINISHED(25); COUNT_UNSET_BIT(0); m += prime * 4 + 1; FALLTHROUGH;
      case 26: CHECK_FINISHED(26); COUNT_UNSET_BIT(6); m += prime * 2 + 1; FALLTHROUGH;
      case 27: CHECK_FINISHED(27); COUNT_UNSET_BIT(5); m += prime * 4 + 2; FALLTHROUGH;
      case 28: CHECK_FINISHED(28); COUNT_UNSET_BIT(2); m += prime * 2 + 1; FALLTHROUGH;
      case 29: CHECK_FINISHED(29); COUNT_UNSET_BIT(1); m += prime * 4 + 1; FALLTHROUGH;
      case 30: CHECK_FINISHED(30); COUNT_UNSET_BIT(7); m += prime * 6 + 3; FALLTHROUGH;
      case 31: CHECK_FINISHED(31); COUNT_UNSET_BIT(4); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 32: CHECK_FINISHED(32); COUNT_UNSET_BIT(4); m += prime * 6 + 3; FALLTHROUGH;
      case 33: CHECK_FINISHED(33); COUNT_UNSET_BIT(7); m += prime * 4 + 3; FALLTHROUGH;
      case 34: CHECK_FINISHED(34); COUNT_UNSET_BIT(1); m += prime * 2 + 1; FALLTHROUGH;
      case 35: CHECK_FINISHED(35); COUNT_UNSET_BIT(2); m += prime * 4 + 2; FALLTHROUGH;
      case 36: CHECK_FINISHED(36); COUNT_UNSET_BIT(5); m += prime * 2 + 1; FALLTHROUGH;
      case 37: CHECK_FINISHED(37); COUNT_UNSET_BIT(6); m += prime * 4 + 3; FALLTHROUGH;
      case 38: CHECK_FINISHED(38); COUNT_UNSET_BIT(0); m += prime * 6 + 3; FALLTHROUGH;
      case 39: CHECK_FINISHED(39); COUNT_UNSET_BIT(3); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 40: CHECK_FINISHED(40); COUNT_UNSET_BIT(5); m += prime * 6 + 4; FALLTHROUGH;
      case 41: CHECK_FINISHED(41); COUNT_UNSET_BIT(3); m += prime * 4 + 2; FALLTHROUGH;
      case 42: CHECK_FINISHED(42); COUNT_UNSET_BIT(7); m += prime * 2 + 2; FALLTHROUGH;
      case 43: CHECK_FINISHED(43); COUNT_UNSET_BIT(1); m += prime * 4 + 2; FALLTHROUGH;
      case 44: CHECK_FINISHED(44); COUNT_UNSET_BIT(6); m += prime * 2 + 2; FALLTHROUGH;
      case 45: CHECK_FINISHED(45); COUNT_UNSET_BIT(0); m += prime * 4 + 2; FALLTHROUGH;
      case 46: CHECK_FINISHED(46); COUNT_UNSET_BIT(4); m += prime * 6 + 4; FALLTHROUGH;
      case 47: CHECK_FINISHED(47); COUNT_UNSET_BIT(2); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 48: CHECK_FINISHED(48); COUNT_UNSET_BIT(6); m += prime * 6 + 5; FALLTHROUGH;
      case 49: CHECK_FINISHED(49); COUNT_UNSET_BIT(2); m += prime * 4 + 3; FALLTHROUGH;
      case 50: CHECK_FINISHED(50); COUNT_UNSET_BIT(3); m += prime * 2 + 1; FALLTHROUGH;
      case 51: CHECK_FINISHED(51); COUNT_UNSET_BIT(7); m += prime * 4 + 4; FALLTHROUGH;
      case 52: CHECK_FINISHED(52); COUNT_UNSET_BIT(0); m += prime * 2 + 1; FALLTHROUGH;
      case 53: CHECK_FINISHED(53); COUNT_UNSET_BIT(4); m += prime * 4 + 3; FALLTHROUGH;
      case 54: CHECK_FINISHED(54); COUNT_UNSET_BIT(5); m += prime * 6 + 5; FALLTHROUGH;
      case 55: CHECK_FINISHED(55); COUNT_UNSET_BIT(1); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 56: CHECK_FINISHED(56); COUNT_UNSET_BIT(7); m += prime * 6 + 6; FALLTHROUGH;
      case 57: CHECK_FINISHED(57); COUNT_UNSET_BIT(6); m += prime * 4 + 4; FALLTHROUGH;
      case 58: CHECK_FINISHED(58); COUNT_UNSET_BIT(5); m += prime * 2 + 2; FALLTHROUGH;
      case 59: CHECK_FINISHED(59); COUNT_UNSET_BIT(4); m += prime * 4 + 4; FALLTHROUGH;
      case 60: CHECK_FINISHED(60); COUNT_UNSET_BIT(3); m += prime * 2 + 2; FALLTHROUGH;
      case 61: CHECK_FINISHED(61); COUNT_UNSET_BIT(2); m += prime * 4 + 4; FALLTHROUGH;
      case 62: CHECK_FINISHED(62); COUNT_UNSET_BIT(1); m += prime * 6 + 6; FALLTHROUGH;
      case 63: CHECK_FINISHED(63); COUNT_UNSET_BIT(0); m += prime * 2 + 1;
    }

    default: UNREACHABLE;
  }
}

} // namespace
