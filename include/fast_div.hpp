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

namespace {

/// Used for (64-bit / 32-bit) = 64-bit.
template <typename X, typename Y>
ALWAYS_INLINE typename std::enable_if<(sizeof(X) == sizeof(uint64_t) &&
                                       sizeof(Y) <= sizeof(uint32_t)), X>::type
fast_div(X x, Y y)
{
  ASSERT(x >= 0);
  ASSERT(y > 0);

  // Unsigned integer division is usually
  // faster than signed integer division.
  using UX = typename pstd::make_unsigned<X>::type;
  using UY = typename pstd::make_unsigned<Y>::type;

#if defined(__x86_64__) && \
   (defined(__GNUC__) || defined(__clang__))

  uint32_t high = uint32_t(UX(x) >> 32);

  if (high < UY(y))
  {
    uint32_t low = uint32_t(x);
    uint32_t d = y;

    // (64-bit / 32-bit) = 32-bit.
    // When we know the result fits into 32-bit (even
    // though the numerator is 64-bit) we can use the divl
    // instruction instead of doing a full 64-bit division.
    __asm__("divl %[divider]"
            : "+a"(low), "+d"(high) : [divider] "r"(d));

    return low;
  }
#endif

  return UX(x) / UY(y);
}

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
  uint64_t high = uint64_t(UX(x) >> 64);

#if defined(__x86_64__) && \
   (defined(__GNUC__) || defined(__clang__))

  if (high < UY(y))
  {
    uint64_t low = uint64_t(x);
    uint64_t d = y;

    // (128-bit / 64-bit) = 64-bit.
    // When we know the result fits into 64-bit (even
    // though the numerator is 128-bit) we can use the divq
    // instruction instead of doing a full 128-bit division.
    __asm__("div %[divider]"
            : "+a"(low), "+d"(high) : [divider] "r"(d));

    return low;
  }
#else
  // This optimization is very important on non x64 CPUs
  // such as ARM64. On AWS Graviton 4 CPUs it e.g. improves
  // performance by about 60% when computing AC(1e22).
  if (high == 0)
    return uint64_t(x) / UY(y);
#endif

  return UX(x) / UY(y);
}

/// Used for  (64-bit /  64-bit) =  64-bit.
/// Used for (128-bit / 128-bit) = 128-bit.
template <typename X, typename Y>
ALWAYS_INLINE typename std::enable_if<(sizeof(X) >= sizeof(uint64_t) &&
                                       sizeof(Y) == sizeof(X)), X>::type
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

#if defined(__x86_64__) && \
   (defined(__GNUC__) || defined(__clang__))

  using UX = typename pstd::make_unsigned<X>::type;

  uint64_t low = uint64_t(x);
  uint64_t high = uint64_t(UX(x) >> 64);
  uint64_t d = y;

  // (128-bit / 64-bit) = 64-bit.
  // When we know the result fits into 64-bit (even
  // though the numerator is 128-bit) we can use the divq
  // instruction instead of doing a full 128-bit division.
  __asm__("div %[divider]"
          : "+a"(low), "+d"(high) : [divider] "r"(d));

  return low;
#else
  return (uint64_t) fast_div(x, y);
#endif
}

/// Used for (64-bit / 32-bit) = 64-bit.
/// Used for (64-bit / 64-bit) = 64-bit.
template <typename X, typename Y>
ALWAYS_INLINE typename std::enable_if<(sizeof(X) <= sizeof(uint64_t) &&
                                       sizeof(Y) <= sizeof(X)), uint64_t>::type
fast_div64(X x, Y y)
{
  return (uint64_t) fast_div(x, y);
}

} // namespace

#endif
