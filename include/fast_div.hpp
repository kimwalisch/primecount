///
/// @file  fast_div.hpp
/// @brief Integer division of small types is much faster than integer
///        division of large types on most CPUs. The fast_div(x, y)
///        function tries to take advantage of this by casting x and y
///        to smaller types (if possible) before doing the division.
///
///        If ENABLE_DIV32 is defined we check at runtime if the
///        dividend and divisor are < 2^32 and if so we use 32-bit
///        integer division instead of 64-bit integer division. On
///        most CPUs before 2020 this significantly improves
///        performance.
///
///        On some new CPUs (such as Intel Cannonlake) 64-bit integer
///        division has been improved significantly and runs as fast
///        as 32-bit integer division. For such CPUs it is best to
///        disable ENABLE_DIV32 (using cmake -DWITH_DIV32=OFF) as this
///        avoids runtime checks for (64-bit / 32-bit) divisions.
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef FAST_DIV_HPP
#define FAST_DIV_HPP

#include <macros.hpp>
#include <int128_t.hpp>

#include <stdint.h>

namespace {

/// Used for (64-bit / 32-bit) = 64-bit.
template <typename X, typename Y>
ALWAYS_INLINE typename std::enable_if<(sizeof(X) == sizeof(uint64_t) &&
                                       sizeof(Y) <= sizeof(uint32_t)), X>::type
fast_div(X x, Y y)
{
  ASSERT(x >= 0);
  ASSERT(y > 0);

#if defined(ENABLE_DIV32)

  using UX = typename pstd::make_unsigned<X>::type;
  using UY = typename pstd::make_unsigned<Y>::type;

  if (x <= pstd::numeric_limits<uint32_t>::max())
    return uint32_t(x) / UY(y);
  else
    return UX(x) / UY(y);
#else
  // Unsigned integer division is usually
  // faster than signed integer division.
  using UX = typename pstd::make_unsigned<X>::type;
  using UY = typename pstd::make_unsigned<Y>::type;
  return UX(x) / UY(y);
#endif
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

  if (x <= pstd::numeric_limits<uint64_t>::max())
    return uint64_t(x) / UY(y);
  else
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

  uint64_t x0 = (uint64_t) x;
  uint64_t x1 = ((uint64_t*) &x)[1];
  uint64_t d = y;

  // (128-bit / 64-bit) = 64-bit.
  // When we know the result fits into 64-bit (even
  // though the numerator is 128-bit) we can use the divq
  // instruction instead of doing a full 128-bit division.
  __asm__("divq %[divider]"
          : "+a"(x0), "+d"(x1) : [divider] "r"(d));

  return x0;
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
