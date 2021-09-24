///
/// @file   isqrt_constexpr.cpp
/// @brief  Test compile time square root function.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <isqrt.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <iostream>

using std::numeric_limits;
using namespace primecount;

#if defined(BAD_ISQRT)

/// The following compile time integer square root function
/// has a recursion depth of O(sqrt(n)). This is very bad, the
/// stack will explose if you try to compute the square root
/// of a number > 10^9. Furthermore constexpr recursion depth
/// is limited by the compiler even more e.g. both GCC and
/// Clang currently limit constexpr recursion depth to 512.
///
/// The nasty thing is that GCC and Clang start calculating
/// constexpr functions at compile time but once the function
/// reaches a recursion depth > 512, GCC and Clang will stop
/// the compile time calculation and instead compute the
/// function at run-time.
///
/// Because of this you cannot easily tell if your constexpr
/// function was evaluated at compile time or at run-time.
/// But you can use static_assert to test whether your constexpr
/// function is calculated at compile time or at runtime.
/// The code below will not compile with the compiler error
/// message telling you the maximum recursion depth has been
/// exceeded.

template <typename T>
constexpr T bad_isqrt_helper(T sq, T dlt, T value)
{
  return sq > value ? (dlt >> 1) - 1
      : bad_isqrt_helper(sq + dlt, dlt + 2, value);
}

template <typename T>
constexpr T bad_isqrt(T value)
{
  return bad_isqrt_helper<T>(1, 3, value);
}

static_assert(bad_isqrt(100000000) == 10000, "bad_isqrt(10^8) failed!");

#endif

int main()
{
  static_assert(ct_sqrt(0) == 0, "ct_sqrt(0) failed!");
  static_assert(ct_sqrt(1) == 1, "ct_sqrt(1) failed!");
  static_assert(ct_sqrt(2) == 1, "ct_sqrt(2) failed!");
  static_assert(ct_sqrt(3) == 1, "ct_sqrt(3) failed!");
  static_assert(ct_sqrt(4) == 2, "ct_sqrt(4) failed!");
  static_assert(ct_sqrt(5) == 2, "ct_sqrt(5) failed!");
  static_assert(ct_sqrt(6) == 2, "ct_sqrt(6) failed!");
  static_assert(ct_sqrt(7) == 2, "ct_sqrt(7) failed!");
  static_assert(ct_sqrt(8) == 2, "ct_sqrt(8) failed!");
  static_assert(ct_sqrt(9) == 3, "ct_sqrt(9) failed!");
  static_assert(ct_sqrt(10) == 3, "ct_sqrt(10) failed!");
  static_assert(ct_sqrt(11) == 3, "ct_sqrt(11) failed!");
  static_assert(ct_sqrt(12) == 3, "ct_sqrt(12) failed!");
  static_assert(ct_sqrt(13) == 3, "ct_sqrt(13) failed!");
  static_assert(ct_sqrt(14) == 3, "ct_sqrt(14) failed!");
  static_assert(ct_sqrt(15) == 3, "ct_sqrt(15) failed!");
  static_assert(ct_sqrt(16) == 4, "ct_sqrt(16) failed!");
  static_assert(ct_sqrt(17) == 4, "ct_sqrt(17) failed!");
  static_assert(ct_sqrt(18) == 4, "ct_sqrt(18) failed!");
  static_assert(ct_sqrt(19) == 4, "ct_sqrt(19) failed!");
  static_assert(ct_sqrt(20) == 4, "ct_sqrt(20) failed!");
  static_assert(ct_sqrt(21) == 4, "ct_sqrt(21) failed!");
  static_assert(ct_sqrt(22) == 4, "ct_sqrt(22) failed!");
  static_assert(ct_sqrt(23) == 4, "ct_sqrt(23) failed!");
  static_assert(ct_sqrt(24) == 4, "ct_sqrt(24) failed!");
  static_assert(ct_sqrt(25) == 5, "ct_sqrt(25) failed!");
  static_assert(ct_sqrt(26) == 5, "ct_sqrt(26) failed!");
  static_assert(ct_sqrt(27) == 5, "ct_sqrt(27) failed!");
  static_assert(ct_sqrt(28) == 5, "ct_sqrt(28) failed!");
  static_assert(ct_sqrt(29) == 5, "ct_sqrt(29) failed!");
  static_assert(ct_sqrt(30) == 5, "ct_sqrt(30) failed!");
  static_assert(ct_sqrt(31) == 5, "ct_sqrt(31) failed!");
  static_assert(ct_sqrt(32) == 5, "ct_sqrt(32) failed!");
  static_assert(ct_sqrt(33) == 5, "ct_sqrt(33) failed!");
  static_assert(ct_sqrt(34) == 5, "ct_sqrt(34) failed!");
  static_assert(ct_sqrt(35) == 5, "ct_sqrt(35) failed!");
  static_assert(ct_sqrt(36) == 6, "ct_sqrt(36) failed!");
  static_assert(ct_sqrt(37) == 6, "ct_sqrt(37) failed!");
  static_assert(ct_sqrt(38) == 6, "ct_sqrt(38) failed!");
  static_assert(ct_sqrt(39) == 6, "ct_sqrt(39) failed!");

  static_assert(ct_sqrt(9223372037000249999ull) == 3037000499ull, "ct_sqrt(3037000500^2-1) failed!");
  static_assert(ct_sqrt(9223372037000250000ull) == 3037000500ull, "ct_sqrt(3037000500^2) failed!");
  static_assert(ct_sqrt(9223372037000250001ull) == 3037000500ull, "ct_sqrt(3037000500^2+1) failed!");

  static_assert(ct_sqrt(numeric_limits<int8_t>::max()) == 11, "ct_sqrt(2^7-1) failed!");
  static_assert(ct_sqrt(numeric_limits<uint8_t>::max()) == 15, "ct_sqrt(2^8-1) failed!");
  static_assert(ct_sqrt(numeric_limits<int16_t>::max()) == 181, "ct_sqrt(2^15-1) failed!");
  static_assert(ct_sqrt(numeric_limits<uint16_t>::max()) == 255, "ct_sqrt(2^16-1) failed!");
  static_assert(ct_sqrt(numeric_limits<int32_t>::max()) == 46340, "ct_sqrt(2^31-1) failed!");
  static_assert(ct_sqrt(numeric_limits<uint32_t>::max()) == 65535, "ct_sqrt(2^32-1) failed!");
  static_assert(ct_sqrt(numeric_limits<int64_t>::max()) == 3037000499ll, "ct_sqrt(2^63-1) failed!");
  static_assert(ct_sqrt(numeric_limits<uint64_t>::max()) == 4294967295ull, "ct_sqrt(2^64-1) failed!");

#if defined(HAVE_INT128_T)

  static_assert(ct_sqrt(numeric_limits<int128_t>::max()) == 13043817825332782212ull, "ct_sqrt(2^127-1) failed!");
  static_assert(ct_sqrt(numeric_limits<uint128_t>::max()) == 18446744073709551615ull, "ct_sqrt(2^128-1) failed!");

  // here std::sqrt((double) 443075998594972078030832658571409090) is 1 too small
  static_assert(ct_sqrt((((int128_t) 24019198012642651) << 64) | 15864680554123835074ull) == 665639541039271553ll, "ct_sqrt(443075998594972078030832658571409090) failed!");
  
  // here std::sqrt((double) 443075998594972075382716071791084150) is 1 too large
  static_assert(ct_sqrt((((int128_t) 24019198012642651) << 64) | 13216563967343510134ull) == 665639541039271551ll, "ct_sqrt(443075998594972075382716071791084150) failed!");
  
  // here std::sqrt((double) 443075998594971958032420320541208365) is 38 too small
  static_assert(ct_sqrt((((int128_t) 24019198012642645) << 64) | 6546732658350944045ull) == 665639541039271462ll, "ct_sqrt(443075998594971958032420320541208365) failed!");

  // here std::sqrt((double) 443075998594971969939937761777907585) is 81 too large
  static_assert(ct_sqrt((((int128_t) 24019198012642646) << 64) | 7506025878091649ull) == 665639541039271471ll, "ct_sqrt(443075998594971969939937761777907585) failed!");

#endif

  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
