///
/// @file   ctz.hpp
/// @brief  Count the number of trailing zeros.
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef CTZ_HPP
#define CTZ_HPP

#include "macros.hpp"
#include <stdint.h>

#if defined(__GNUC__) || \
    __has_builtin(__builtin_ctzl)

namespace {

inline int ctz64(uint64_t x)
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

} // namespace

#elif __cplusplus >= 202002L && \
      __has_include(<bit>)

#include <bit>

namespace {

inline int ctz64(uint64_t x)
{
  return std::countr_zero(x);
}

} // namespace

#elif defined(_MSC_VER) && \
      __has_include(<intrin.h>)

#include <intrin.h>

namespace {

inline int ctz64(uint64_t x)
{
  // _BitScanForward64(0) is undefined behavior
  ASSERT(x != 0);

  unsigned long r;
  _BitScanForward64(&r, x);
  return (int) r;
}

} // namespace

#endif

#endif // CTZ_HPP
