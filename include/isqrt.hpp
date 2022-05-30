///
/// @file  isqrt.hpp
/// @brief Integer square root function
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef ISQRT_HPP
#define ISQRT_HPP

#include <macros.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdint.h>
#include <type_traits>

namespace {

#if __cplusplus >= 202002L

/// C++20 compile time square root using binary search
template <typename T>
consteval T sqrt_helper(T x, T lo, T hi)
{
  if (lo == hi)
    return lo;

  const T mid = (lo + hi + 1) / 2;

  if (x / mid < mid)
    return sqrt_helper<T>(x, lo, mid - 1);
  else
    return sqrt_helper(x, mid, hi);
}

template <typename T>
consteval T ct_sqrt(T x)
{
  return sqrt_helper<T>(x, 0, x / 2 + 1);
}

#elif __cplusplus >= 201402L

/// C++14 compile time square root using binary search
template <typename T>
constexpr T sqrt_helper(T x, T lo, T hi)
{
  if (lo == hi)
    return lo;

  const T mid = (lo + hi + 1) / 2;

  if (x / mid < mid)
    return sqrt_helper<T>(x, lo, mid - 1);
  else
    return sqrt_helper(x, mid, hi);
}

template <typename T>
constexpr T ct_sqrt(T x)
{
  return sqrt_helper<T>(x, 0, x / 2 + 1);
}

#else

#define MID ((lo + hi + 1) / 2)

/// C++11 compile time square root using binary search
template <typename T>
constexpr T sqrt_helper(T x, T lo, T hi)
{
  return lo == hi ? lo : ((x / MID < MID)
      ? sqrt_helper<T>(x, lo, MID - 1) : sqrt_helper<T>(x, MID, hi));
}

template <typename T>
constexpr T ct_sqrt(T x)
{
  return sqrt_helper<T>(x, 0, x / 2 + 1);
}

#endif

template <typename T>
ALWAYS_INLINE T isqrt(T x)
{
  T s = (T) std::sqrt((double) x);

  // By using constexpr for the sqrt_max variable type it
  // is guaranteed that ct_sqrt() is evaluated at compile
  // time. Compilation will fail if the compiler fails to
  // evaluate ct_sqrt() at compile time. This is great,
  // ct_sqrt() must be evaluated at compile time otherwise
  // the runtime complexity of isqrt(x) would deteriorate
  // from O(1) to O(log2(x)).
  //
  // If sqrt_max were declared without constexpr then the
  // compiler would be free to compute ct_sqrt() either at
  // compile time or at run time e.g. GCC-11 computes
  // ct_sqrt(MAX_INT128) at compile time whereas Clang-12
  // computes ct_sqrt(MAX_INT128) at run time even at -O2.
  //
  // C++20 fixed this annoying issue by adding consteval
  // to C++. Hence if the compiler supports C++20 ct_sqrt()
  // is defined as consteval instead of constexpr. Hence
  // using C++20 ct_sqrt() will be evaluated at compile
  // time in all cases i.e. even if sqrt_max were declared
  // without constexpr.
  //
  constexpr T sqrt_max = ct_sqrt(std::numeric_limits<T>::max());

  // For 128-bit integers we use uint64_t as the
  // result type. For all other types we use the
  // same result type as the input type.
  using R = typename std::conditional<sizeof(T) / 2 == sizeof(uint64_t), uint64_t, T>::type;
  R r = (R) std::min(s, sqrt_max);

  // In my tests the first corrections were needed above
  // 10^22 where the results were off by 1. Above 10^32 the
  // first results occurred that were off by > 1. Since
  // primecount only supports numbers up to 10^31 this is
  // not an issue for us.
  if (r * (T) r > x)
  {
    do { r--; }
    while (r * (T) r > x);
  }
  // Same as (r + 1)^2 < x but overflow safe
  else if ((T) (r * 2) < x - r * (T) r)
  {
    do { r++; }
    while ((T) (r * 2) < x - r * (T) r);
  }

  return r;
}

} // namespace

#endif
