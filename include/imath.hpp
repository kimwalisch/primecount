///
/// @file  imath.hpp
/// @brief Integer math functions
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef IMATH_HPP
#define IMATH_HPP

#include <isqrt.hpp>

#include <stdint.h>
#include <cmath>
#include <limits>

#if __cplusplus >= 202002L
  #include <bit>
  #include <type_traits>
#endif

namespace {

inline int64_t isquare(int64_t x)
{
  return x * x;
}

template <typename A, typename B>
inline A ceil_div(A a, B b)
{
  return (A) ((a + b - 1) / b);
}

/// Next power of 2 >= x
template <typename T>
inline T next_power_of_2(T x)
{
#if __cplusplus >= 202002L
  auto ux = std::make_unsigned_t<T>(x);
  return std::bit_ceil(ux);
#else
  if (x == 0)
    return 1;

  x--;
  T bits = std::numeric_limits<T>::digits;

  for (T i = 1; i < bits; i += i)
    x |= (x >> i);

  return ++x;
#endif
}

template <typename T>
inline int ilog(T x)
{
  return (int) std::log((double) x);
}

template <typename T>
inline T ilog2(T x)
{
#if __cplusplus >= 202002L
  auto ux = std::make_unsigned_t<T>(x);
  ux = (ux > 0) ? ux : 1;
  return std::bit_width(ux) - 1;
#else
  T log2 = 0;
  T bits = std::numeric_limits<T>::digits;

  for (T i = bits / 2; i > 0; i /= 2)
  {
    T one = 1;
    if (x >= (one << i))
    {
      x >>= i;
      log2 += i;
    }
  }

  return log2;
#endif
}

template <typename T>
inline T ipow(T x, int n)
{
  T r = 1;
  for (int i = 0; i < n; i++)
    r *= x;

  return r;
}

/// Integer nth root
template <int N, typename T>
inline T iroot(T x)
{
  T r;

  if (N == 3)
    r = (T) std::cbrt((double) x);
  else
    r = (T) std::pow((double) x, 1.0 / N);

  // fix root too large
  for (; r > 0; r--)
    if (ipow(r, N - 1) <= x / r)
      break;

  // fix root too small
  while (ipow(r + 1, N - 1) <= x / (r + 1))
    r += 1;

  return r;
}

} // namespace

#endif
