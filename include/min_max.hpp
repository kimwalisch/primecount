///
/// @file  min_max.hpp
/// @brief Template min and max functions that allow comparing
///        different types if both types are integral
///        and sizeof(A) >= sizeof(B).
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef MIN_MAX_HPP
#define MIN_MAX_HPP

#include <int128_t.hpp>

#include <algorithm>
#include <type_traits>

namespace primecount {

/// Check if A and B are:
/// 1) Same types
/// 2) Different integral types which satisfy:
///      sizeof(A) > sizeof(B) ||
///     (sizeof(A) == sizeof(B) && is_unsigned<A> && is_signed<B>)
///
template <typename A, typename B>
struct is_comparable
{
  enum {
    value = std::is_same<A, B>::value || ((
            prt::is_integral<A>::value &&
            prt::is_integral<B>::value) && (
            sizeof(A) > sizeof(B) || (
            sizeof(A) == sizeof(B) &&
            std::is_unsigned<A>::value &&
            std::is_signed<B>::value)))
  };
};

template <typename A, typename B>
inline B min(A a, B b)
{
  static_assert(is_comparable<A, B>::value,
                "min(A, B): Cannot compare types A and B");

  return (B) std::min(a, (A) b);
}

template <typename A, typename B>
inline B min3(A a, B b, B c)
{
  static_assert(is_comparable<A, B>::value,
                "min3(A, B, B): Cannot compare types A and B");

  return (B) std::min(a, (A) std::min(b, c));
}

template <typename A, typename B>
inline A max(A a, B b)
{
  static_assert(is_comparable<A, B>::value,
                "max(A, B): Cannot compare types A and B");

  return std::max(a, (A) b);
}

template <typename A, typename B>
inline A max3(A a, B b, B c)
{
  static_assert(is_comparable<A, B>::value,
                "max3(A, B, B): Cannot compare types A and B");

  return std::max(a, (A) std::max(b, c));
}

} // namespace

#endif
