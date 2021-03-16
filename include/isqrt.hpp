///
/// @file  isqrt.hpp
/// @brief Integer square root function
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef ISQRT_HPP
#define ISQRT_HPP

#include <algorithm>
#include <cmath>
#include <limits>

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
inline T isqrt(T x)
{
  T r = (T) std::sqrt((double) x);

  constexpr T sqrt_max = ct_sqrt(std::numeric_limits<T>::max());
  r = std::min(r, sqrt_max);

  while (r * r > x)
    r--;
  while (x - r * r > r * 2)
    r++;

  return r;
}

} // namespace

#endif
