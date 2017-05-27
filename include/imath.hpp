///
/// @file  imath.hpp
/// @brief Integer math functions used in primecount
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef IMATH_HPP
#define IMATH_HPP

#include <isqrt.hpp>

#include <stdint.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

namespace {

inline int64_t isquare(int64_t x)
{
  return x * x;
}

template <typename A, typename B>
inline A ceil_div(A a, B b)
{
  assert(b > 0);
  return (A) ((a + b - 1) / b);
}

template <typename T>
inline T number_of_bits(T)
{
  return (T) (sizeof(T) * 8);
}

template <typename T>
inline bool is_power_of_2(T x)
{
  return (x != 0 && (x & (x - 1)) == 0);
}

template <typename T>
inline T next_power_of_2(T x)
{
  if (x == 0)
    return 1;

  x--;
  for (T i = 1; i < number_of_bits(x); i += i)
    x |= (x >> i);

  return ++x;
}

template <typename T>
inline T prev_power_of_2(T x)
{
  for (T i = 1; i < number_of_bits(x); i += i)
    x |= (x >> i);

  return x - (x >> 1);
}

template <typename T>
inline int ilog(T x)
{
  return (int) std::log((double) x);
}

template <typename T>
inline T ipow(T x, int n)
{
  T r = 1;
  for (int i = 0; i < n; i++)
    r *= x;

  return r;
}

/// Check if ipow(x, n) <= limit
template <typename T>
inline bool ipow_less_equal(T x, int n, T limit)
{
  if (limit <= 0)
    return false;

  for (T r = 1; n > 0; n--, r *= x)
    if (r > limit / x)
      return false;

  return true;
}

/// Integer nth root
template <int N, typename T>
inline T iroot(T x)
{
  T r = (T) std::pow((double) x, 1.0 / N);

  while (ipow(r, N) > x)
    r--;
  while (ipow_less_equal(r + 1, N, x))
    r++;

  return r;
}

/// Count the number of primes <= x using binary search.
/// @pre primes[1] = 2, primes[3] = 3, ...
/// @pre x <= primes.back()
///
template <typename T1, typename T2>
inline T2 pi_bsearch(const std::vector<T1>& primes, T2 x)
{
  assert(primes[0] == 0);
  return (T2) (std::upper_bound(primes.begin() + 1, primes.end(), x) - (primes.begin() + 1));
}

template <typename T1, typename T2, typename T3>
inline T2 in_between(T1 min, T2 x, T3 max)
{
  if (x < min)
    return (T2) min;
  if (x > max)
    return (T2) max;

  return x;
}

} // namespace

#endif
