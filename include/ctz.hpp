///
/// @file   ctz.hpp
/// @brief  Count the number of leading and trailing zeros.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef CTZ_HPP
#define CTZ_HPP

#include "macros.hpp"

#include <climits>
#include <stdint.h>

#if defined(__GNUC__) || \
    __has_builtin(__builtin_ctzl)

namespace primecount {

ALWAYS_INLINE int ctz64(uint64_t x)
{
  // __builtin_ctz(0) is undefined behavior
  ASSERT(x != 0);

#if __cplusplus >= 201703L
  if constexpr(sizeof(int) >= sizeof(uint64_t))
    return __builtin_ctz(x);
  else if constexpr(sizeof(long) >= sizeof(uint64_t))
    return __builtin_ctzl(x);
  else if constexpr(sizeof(long long) >= sizeof(uint64_t))
    return __builtin_ctzll(x);
#else
    return __builtin_ctzll(x);
#endif
}

ALWAYS_INLINE int clz64(uint64_t x)
{
  // __builtin_clz(0) is undefined behavior
  ASSERT(x != 0);

#if __cplusplus >= 201703L
  if constexpr(sizeof(int) >= sizeof(uint64_t))
    return __builtin_clz(x) - int(sizeof(int) - sizeof(uint64_t)) * CHAR_BIT;
  else if constexpr(sizeof(long) >= sizeof(uint64_t))
    return __builtin_clzl(x) - int(sizeof(long) - sizeof(uint64_t)) * CHAR_BIT;
  else if constexpr(sizeof(long long) >= sizeof(uint64_t))
    return __builtin_clzll(x) - int(sizeof(long long) - sizeof(uint64_t)) * CHAR_BIT;
#else
  return __builtin_clzll(x) - int(sizeof(long long) - sizeof(uint64_t)) * CHAR_BIT;
#endif
}

} // namespace

#elif __cplusplus >= 202002L && \
      __has_include(<bit>)

#include <bit>

namespace primecount {

ALWAYS_INLINE int ctz64(uint64_t x)
{
  return std::countr_zero(x);
}

ALWAYS_INLINE int clz64(uint64_t x)
{
  return std::countl_zero(x);
}

} // namespace

#elif defined(_MSC_VER) && \
      defined(_WIN64) && \
      __has_include(<intrin.h>)

#include <intrin.h>

namespace primecount {

ALWAYS_INLINE int ctz64(uint64_t x)
{
  // _BitScanForward64(0) is undefined behavior
  ASSERT(x != 0);

  unsigned long r;
  _BitScanForward64(&r, x);
  return (int) r;
}

ALWAYS_INLINE int clz64(uint64_t x)
{
  // _BitScanReverse64(0) is undefined behavior
  ASSERT(x != 0);

  unsigned long r;
  _BitScanReverse64(&r, x);
  return 63 - (int) r;
}

} // namespace

#elif defined(_MSC_VER) && \
      (defined(_M_IX86) || defined(_M_ARM)) && \
      __has_include(<intrin.h>)

#include <intrin.h>

namespace primecount {

ALWAYS_INLINE int ctz64(uint64_t x)
{
  // _BitScanForward(0) is undefined behavior
  ASSERT(x != 0);

  unsigned long r;
  if (_BitScanForward(&r, (unsigned long) x))
    return (int) r;

  _BitScanForward(&r, (unsigned long) (x >> 32));
  return (int) r + 32;
}

ALWAYS_INLINE int clz64(uint64_t x)
{
  // _BitScanReverse(0) is undefined behavior
  ASSERT(x != 0);

  unsigned long r;
  if (_BitScanReverse(&r, (unsigned long) (x >> 32)))
    return 31 - (int) r;

  _BitScanReverse(&r, (unsigned long) x);
  return 63 - (int) r;
}

} // namespace

#endif

#endif // CTZ_HPP
