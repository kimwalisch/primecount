///
/// @file  fast_div.hpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef FAST_DIV_HPP
#define FAST_DIV_HPP

#include <cassert>
#include <limits>

namespace primecount {

inline int64_t fast_div(int64_t x, int32_t y)
{
  assert(x >= 0);
  assert(y >= 0);

  if (x <= std::numeric_limits<uint32_t>::max())
    return (uint32_t) x / (uint32_t) y;

  return x / y;
}

} // namespace

#endif
