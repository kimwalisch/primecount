///
/// @file  isqrt.hpp
/// @brief Integer square root function
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef ISQRT_HPP
#define ISQRT_HPP

#include <int128_t.hpp>

#include <algorithm>
#include <cmath>

namespace primecount {

template <typename T>
struct isqrt_limits { };

template <> struct isqrt_limits<int8_t>
{ static constexpr int8_t max() { return 11; } };

template <> struct isqrt_limits<uint8_t>
{ static constexpr uint8_t max() { return 15; } };

template <> struct isqrt_limits<int16_t>
{ static constexpr int16_t max() { return 181; } };

template <> struct isqrt_limits<uint16_t>
{ static constexpr uint16_t max() { return 255; } };

template <> struct isqrt_limits<int32_t>
{ static constexpr int32_t max() { return 46340; } };

template <> struct isqrt_limits<uint32_t>
{ static constexpr uint32_t max() { return 65535; } };

template <> struct isqrt_limits<int64_t>
{ static constexpr int64_t max() { return 3037000499ll; } };

template <> struct isqrt_limits<uint64_t>
{ static constexpr uint64_t max() { return 4294967295ull; } };

template <> struct isqrt_limits<int128_t>
{ static constexpr uint64_t max() { return 13043817825332782212ull; } };

template <> struct isqrt_limits<uint128_t>
{ static constexpr uint64_t max() { return 18446744073709551615ull; } };

template <typename T>
inline T isqrt(T x)
{
  T r = (T) std::sqrt((double) x);

  constexpr T max_sqrt = isqrt_limits<T>::max();
  r = std::min(r, max_sqrt);

  while (r * r > x)
    r--;
  while (x - r * r > r * 2)
    r++;

  return r;
}

} // namespace

#endif
