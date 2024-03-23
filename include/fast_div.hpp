///
/// @file  fast_div.hpp
/// @brief Integer division of small types is much faster than integer
///        division of large types on most CPUs. The fast_div(x, y)
///        function tries to take advantage of this by casting x and y
///        to smaller types (if possible) before doing the division.
///
///        If ENABLE_DOUBLE_INTEGER_DIVISION is defined, we check at
///        runtime if the dividend and divisor are <= 2^53 (i.e. these
///        integers can be represented exactly by the double type)
///        and if so we use double floating-point division instead of
///        integer division and convert the result back to an integer.
///        This significantly improves performance on many CPU
///        architectures where floating-point division usually runs
///        much faster than integer division.
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef FAST_DIV_HPP
#define FAST_DIV_HPP

#include <macros.hpp>

#include <limits>
#include <stdint.h>
#include <type_traits>

#if defined(ENABLE_DOUBLE_INTEGER_DIVISION)

static_assert(std::numeric_limits<double>::radix == 2,
  "double type radix != 2");
static_assert(std::numeric_limits<double>::digits <= std::numeric_limits<unsigned long long>::digits,
  "double mantissa bits > long long bits!");

#endif

namespace {

/// Converts an 128-bit integer type to uint64_t, other smaller
/// integer types are simply converted to unsigned.
/// 
/// make_smaller<uint128_t>::type -> uint64_t.
/// make_smaller< uint64_t>::type -> uint64_t.
/// make_smaller< uint32_t>::type -> uint32_t.
///
template <typename T>
struct make_smaller
{
  using type = typename std::make_unsigned<
                 typename std::conditional<
                   sizeof(T) == sizeof(uint64_t) * 2,
                     uint64_t, T>::type>::type;
};

/// Used for (64-bit / 32-bit) = 64-bit.
template <typename X, typename Y>
ALWAYS_INLINE constexpr
typename std::enable_if<(sizeof(X) == sizeof(uint64_t) &&
                         sizeof(Y) <= sizeof(uint32_t)), X>::type
fast_div(X x, Y y)
{
  ASSERT(x >= 0);
  ASSERT(y > 0);

#if defined(ENABLE_DOUBLE_INTEGER_DIVISION)

  using UX = typename std::make_unsigned<X>::type;
  using UY = typename std::make_unsigned<Y>::type;

  if (UX(x) <= (1ull << std::numeric_limits<double>::digits))
    return X((double) x / (double) y);
  else
    return UX(x) / UY(y);
#else
  // Unsigned integer division is usually
  // faster than signed integer division.
  using UX = typename std::make_unsigned<X>::type;
  using UY = typename std::make_unsigned<Y>::type;
  return UX(x) / UY(y);
#endif
}

/// Used for  (64-bit /  64-bit) =  64-bit.
/// Used for (128-bit / 128-bit) = 128-bit.
template <typename X, typename Y>
ALWAYS_INLINE constexpr
typename std::enable_if<(sizeof(X) >= sizeof(uint64_t) &&
                         sizeof(Y) == sizeof(X)), X>::type
fast_div(X x, Y y)
{
  ASSERT(x >= 0);
  ASSERT(y > 0);

  // Unsigned integer division is usually
  // faster than signed integer division.
  using UX = typename std::make_unsigned<X>::type;
  using UY = typename std::make_unsigned<Y>::type;
  return UX(x) / UY(y);
}

/// Used for (128-bit / 32-bit) = 128-bit.
/// Used for (128-bit / 64-bit) = 128-bit.
template <typename X, typename Y>
ALWAYS_INLINE constexpr
typename std::enable_if<(sizeof(X) > sizeof(uint64_t) &&
                         sizeof(Y) <= sizeof(uint64_t)), X>::type
fast_div(X x, Y y)
{
  ASSERT(x >= 0);
  ASSERT(y > 0);

  // Unsigned integer division is usually
  // faster than signed integer division.
  using UX = typename std::make_unsigned<X>::type;
  using UY = typename std::make_unsigned<Y>::type;
  using smaller_t = typename make_smaller<X>::type;

  if (x <= std::numeric_limits<smaller_t>::max())
    return smaller_t(x) / UY(y);
  else
    return UX(x) / UY(y);
}

/// Used for (128-bit / 32-bit) = 64-bit.
/// Used for (128-bit / 64-bit) = 64-bit.
/// Use this function only when you know for sure
/// that the result is < 2^64.
///
template <typename X, typename Y>
ALWAYS_INLINE
typename std::enable_if<(sizeof(X) > sizeof(uint64_t) &&
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
ALWAYS_INLINE
typename std::enable_if<(sizeof(X) <= sizeof(uint64_t) &&
                         sizeof(Y) <= sizeof(X)), uint64_t>::type
fast_div64(X x, Y y)
{
  return (uint64_t) fast_div(x, y);
}

} // namespace

#endif
