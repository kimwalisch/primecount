///
/// @file  min.hpp
/// @brief Template min and max functions that allow comparing
///        different types if both types are integral
///        and sizeof(A) >= sizeof(B).
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef MIN_HPP
#define MIN_HPP

#include <int128_t.hpp>
#include <macros.hpp>

#include <algorithm>

namespace {

template <typename A, typename B>
struct is_comparable
{
  enum {
    value = pstd::is_same<A, B>::value || (
            pstd::is_integral<A>::value &&
            pstd::is_integral<B>::value &&
            sizeof(A) >= sizeof(B))
  };
};

template <typename A, typename B>
ALWAYS_INLINE B min(A a, B b)
{
  static_assert(is_comparable<A, B>::value,
                "min(A, B): Cannot compare types A and B");

#if defined(ENABLE_ASSERT)
  if (pstd::is_unsigned<A>::value && pstd::is_signed<B>::value) ASSERT(b >= 0);
  if (pstd::is_unsigned<B>::value && pstd::is_signed<A>::value) ASSERT(a >= 0 && b <= pstd::numeric_limits<A>::max());
#endif

  return (B) std::min(a, (A) b);
}

template <typename A, typename B>
ALWAYS_INLINE A max(A a, B b)
{
  static_assert(is_comparable<A, B>::value,
                "max(A, B): Cannot compare types A and B");

#if defined(ENABLE_ASSERT)
  if (pstd::is_unsigned<A>::value && pstd::is_signed<B>::value) ASSERT(b >= 0);
  if (pstd::is_unsigned<B>::value && pstd::is_signed<A>::value) ASSERT(b <= pstd::numeric_limits<A>::max());
#endif

  return std::max(a, (A) b);
}

template <typename A, typename B, typename C>
ALWAYS_INLINE C min3(A a, B b, C c)
{
  return min(a, min(b, c));
}

template <typename A, typename B, typename C>
ALWAYS_INLINE A max3(A a, B b, C c)
{
  return max(a, max(b, c));
}

} // namespace

#endif
