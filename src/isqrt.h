///
/// @file  isqrt.h
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef ISQRT_PRIMECOUNT_H
#define ISQRT_PRIMECOUNT_H

#include <stdint.h>
#include <cmath>

namespace primecount {

inline int32_t isqrt(int64_t x)
{
  double x2 = static_cast<double>(x);
  int64_t r = static_cast<int64_t>(std::sqrt(x2));
  // correct rounding errors
  while (r * r > x)
    r--;
  while ((r + 1) * (r + 1) <= x)
    r++;
  return static_cast<int32_t>(r);
}

inline int32_t isqrt3(int64_t x)
{
  double x2 = static_cast<double>(x);
  int64_t r = static_cast<int64_t>(std::pow(x2, 1.0 / 3.0));
  // correct rounding errors
  while (r * r * r > x)
    r--;
  while ((r + 1) * (r + 1) * (r + 1) <= x)
    r++;
  return static_cast<int32_t>(r);
}

inline int32_t isqrt4(int64_t x)
{
  double x2 = static_cast<double>(x);
  int64_t r = static_cast<int64_t>(std::pow(x2, 1.0 / 4.0));
  // correct rounding errors
  while (r * r * r * r > x)
    r--;
  while ((r + 1) * (r + 1) * (r + 1) * (r + 1) <= x)
    r++;
  return static_cast<int32_t>(r);
}

} // namespace primecount

#endif
