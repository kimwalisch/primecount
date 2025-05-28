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
      defined(_WIN64) && \
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

#else

// Portable pure integer count trailing zeros algorithm.
// https://www.chessprogramming.org/BitScan#With_separated_LS1B

namespace {

const int index64[64] =
{
  0, 47,  1, 56, 48, 27,  2, 60,
  57, 49, 41, 37, 28, 16,  3, 61,
  54, 58, 35, 52, 50, 42, 21, 44,
  38, 32, 29, 23, 17, 11,  4, 62,
  46, 55, 26, 59, 40, 36, 15, 53,
  34, 51, 20, 43, 31, 22, 10, 45,
  25, 39, 14, 33, 19, 30,  9, 24,
  13, 18,  8, 12,  7,  6,  5, 63
};

inline int ctz64(uint64_t x)
{
  // ctz64(0) is undefined behavior
  ASSERT(x != 0);
  constexpr uint64_t debruijn64 = 0x03f79d71b4cb0a89ull;
  return index64[((x ^ (x-1)) * debruijn64) >> 58];
}

} // namespace

#endif

#endif // CTZ_HPP
