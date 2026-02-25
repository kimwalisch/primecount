///
/// @file  fast_div.hpp
/// @brief Integer division of small types is much faster than integer
///        division of large types on most CPUs. The fast_div(x, y)
///        function tries to take advantage of this by casting x and y
///        to smaller types (if possible) before doing the division.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef FAST_DIV_HPP
#define FAST_DIV_HPP

#include <macros.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <type_traits>

#if defined(HAVE_INT128_T)
  #include <fast_div_128_64.hpp>
#endif

namespace {

/// Used for (64-bit / 32-bit) = 64-bit.
/// Used for (64-bit / 64-bit) = 64-bit.
template <typename X, typename Y>
ALWAYS_INLINE typename std::enable_if<(sizeof(X) <= sizeof(uint64_t) &&
                                       sizeof(Y) <= sizeof(uint64_t)), X>::type
fast_div(X x, Y y)
{
  ASSERT(x >= 0);
  ASSERT(y > 0);

  // Unsigned integer division is usually
  // faster than signed integer division.
  using UX = typename pstd::make_unsigned<X>::type;
  using UY = typename pstd::make_unsigned<Y>::type;
  return UX(x) / UY(y);
}

#if defined(HAVE_INT128_T)

/// Used for (128-bit / 32-bit) = 128-bit.
/// Used for (128-bit / 64-bit) = 128-bit.
template <typename X, typename Y>
ALWAYS_INLINE typename std::enable_if<(sizeof(X) > sizeof(uint64_t) &&
                                       sizeof(Y) <= sizeof(uint64_t)), X>::type
fast_div(X x, Y y)
{
  ASSERT(x >= 0);
  ASSERT(y > 0);

  // Unsigned integer division is usually
  // faster than signed integer division.
  using UX = typename pstd::make_unsigned<X>::type;
  using UY = typename pstd::make_unsigned<Y>::type;
  uint64_t numhi = uint64_t(UX(x) >> 64);

  if (numhi < uint64_t(y))
    return fast_div_128_to_64(x, y);
  else
    return UX(x) / UY(y);
}

#endif

/// Used for (64-bit / 32-bit) = 64-bit.
/// Used for (64-bit / 64-bit) = 64-bit.
template <typename X, typename Y>
ALWAYS_INLINE typename std::enable_if<(sizeof(X) <= sizeof(uint64_t) &&
                                       sizeof(Y) <= sizeof(X)), uint64_t>::type
fast_div64(X x, Y y)
{
  ASSERT(x >= 0);
  ASSERT(y > 0);

  // Unsigned integer division is usually
  // faster than signed integer division.
  using UX = typename pstd::make_unsigned<X>::type;
  using UY = typename pstd::make_unsigned<Y>::type;
  return UX(x) / UY(y);
}

#if defined(HAVE_INT128_T)

/// Used for (128-bit / 32-bit) = 64-bit.
/// Used for (128-bit / 64-bit) = 64-bit.
/// Use this function only when you know for sure
/// that the result is < 2^64.
///
template <typename X, typename Y>
ALWAYS_INLINE typename std::enable_if<(sizeof(X) > sizeof(uint64_t) &&
                                       sizeof(Y) <= sizeof(uint64_t)), uint64_t>::type
fast_div64(X x, Y y)
{
  ASSERT(x >= 0);
  ASSERT(y > 0);

  return fast_div_128_to_64(x, y);
}

#endif

} // namespace

#endif
