///
/// @file   int128_STL_patch.hpp
/// @brief  Add missing int128_t support to the C++ STL.
///         This is required using GCC/Clang if the user compiles with
///         -std=c++11 instead of -std=gnu++11 and also for LLVM/Clang
///         on Windows.
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef INT128_STL_PATCH_HPP
#define INT128_STL_PATCH_HPP

#include <limits>
#include <type_traits>

namespace std {

/// std::is_integral
template<> struct is_integral<primecount::int128_t> : std::true_type {};
template<> struct is_integral<primecount::uint128_t> : std::true_type {};

/// std::is_floating_point
template<> struct is_floating_point<primecount::int128_t> : std::false_type {};
template<> struct is_floating_point<primecount::uint128_t> : std::false_type {};

/// std::is_signed
template<> struct is_signed<primecount::int128_t> : std::true_type {};
template<> struct is_signed<primecount::uint128_t> : std::false_type {};

/// std::is_unsigned
template<> struct is_unsigned<primecount::int128_t> : std::false_type {};
template<> struct is_unsigned<primecount::uint128_t> : std::true_type {};

/// std::make_unsigned
template<> struct make_unsigned<primecount::int128_t> { using type = primecount::uint128_t; };
template<> struct make_unsigned<primecount::uint128_t> { using type = primecount::uint128_t; };

/// std::numeric_limits
template<> struct numeric_limits<primecount::int128_t>
{
  static constexpr primecount::int128_t min() { return primecount::int128_t(primecount::uint128_t(1) << 127); }
  static constexpr primecount::int128_t max() { return primecount::int128_t((primecount::uint128_t(1) << 127) - 1); }
  static constexpr int digits = 127;
};

/// std::numeric_limits
template<> struct numeric_limits<primecount::uint128_t>
{
  static constexpr primecount::uint128_t min() { return 0; }
  static constexpr primecount::uint128_t max() { return ~primecount::uint128_t(0); }
  static constexpr int digits = 128;
};

} // namespace

#endif
