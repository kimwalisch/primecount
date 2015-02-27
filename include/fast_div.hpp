///
/// @file  fast_div.hpp
/// @brief Integer division of small types is much faster than
///        integer division of large types on most CPUs. The
///        fast_div(x, y) function tries to take advantage of this
///        by casting x and y to smaller types (if possible) before
///        doing the division.
///
/// Copyright (C) 2015 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef FAST_DIV_HPP
#define FAST_DIV_HPP

#include <int128.hpp>
#include <cassert>
#include <limits>

namespace primecount {

// x = 64-bit and y = 32-bit

inline int64_t fast_div(int64_t x, int32_t y)
{
  assert(x >= 0);
  assert(y >= 0);

  if (x <= std::numeric_limits<uint32_t>::max())
    return (uint32_t) x / (uint32_t) y;

  return x / y;
}

inline int64_t fast_div(int64_t x, uint32_t y)
{
  assert(x >= 0);

  if (x <= std::numeric_limits<uint32_t>::max())
    return (uint32_t) x /  y;

  return x / y;
}

inline uint64_t fast_div(uint64_t x, int32_t y)
{
  assert(y >= 0);

  if (x <= std::numeric_limits<uint32_t>::max())
    return (uint32_t) x /  y;

  return x / y;
}

inline uint64_t fast_div(uint64_t x, uint32_t y)
{
  if (x <= std::numeric_limits<uint32_t>::max())
    return (uint32_t) x /  y;

  return x / y;
}

// x = 64-bit and y = 64-bit

inline int64_t fast_div(int64_t x, int64_t y)
{
  return x / y;
}

inline uint64_t fast_div(int64_t x, uint64_t y)
{
  return x / y;
}

inline uint64_t fast_div(uint64_t x, int64_t y)
{
  return x / y;
}

inline uint64_t fast_div(uint64_t x, uint64_t y)
{
  return x / y;
}

#if defined(HAVE_INT128_T)

// x = 128-bit and y = 32-bit

inline int128_t fast_div(int128_t x, int32_t y)
{
  assert(x >= 0);
  assert(y >= 0);

  if (x <= std::numeric_limits<uint64_t>::max())
    return (uint64_t) x / (uint64_t) y;

  return x / y;
}

inline int128_t fast_div(int128_t x, uint32_t y)
{
  assert(x >= 0);

  if (x <= std::numeric_limits<uint64_t>::max())
    return (uint64_t) x / (uint64_t) y;

  return x / y;
}

inline uint128_t fast_div(uint128_t x, int32_t y)
{
  assert(y >= 0);

  if (x <= std::numeric_limits<uint64_t>::max())
    return (uint64_t) x / (uint64_t) y;

  return x / y;
}

inline uint128_t fast_div(uint128_t x, uint32_t y)
{
  if (x <= std::numeric_limits<uint64_t>::max())
    return (uint64_t) x / (uint64_t) y;

  return x / y;
}

// x = 128-bit and y = 64-bit

inline int128_t fast_div(int128_t x, int64_t y)
{
  assert(x >= 0);
  assert(y >= 0);

  if (x <= std::numeric_limits<uint64_t>::max())
    return (uint64_t) x / (uint64_t) y;

  return x / y;
}

inline int128_t fast_div(int128_t x, uint64_t y)
{
  assert(x >= 0);

  if (x <= std::numeric_limits<uint64_t>::max())
    return (uint64_t) x / (uint64_t) y;

  return x / y;
}

inline uint128_t fast_div(uint128_t x, int64_t y)
{
  assert(y >= 0);

  if (x <= std::numeric_limits<uint64_t>::max())
    return (uint64_t) x / (uint64_t) y;

  return x / y;
}

inline uint128_t fast_div(uint128_t x, uint64_t y)
{
  if (x <= std::numeric_limits<uint64_t>::max())
    return (uint64_t) x / (uint64_t) y;

  return x / y;
}

// x = 128-bit and y = 128-bit

inline int128_t fast_div(int128_t x, int128_t y)
{
  return x / y;
}

inline uint128_t fast_div(int128_t x, uint128_t y)
{
  return x / y;
}

inline uint128_t fast_div(uint128_t x, int128_t y)
{
  return x / y;
}

inline uint128_t fast_div(uint128_t x, uint128_t y)
{
  return x / y;
}

#endif /* HAVE_INT128_T */

} // namespace

#endif
