///
/// @file   int128.hpp
/// @brief  Convenience functions and typedefs for the non standard
///         __int128_t & __uint128_t integer types.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PTYPES_HPP
#define PTYPES_HPP

#include <stdint.h>

#if defined(HAVE_INT128_T)

namespace primecount {

typedef int128_t maxint_t;
typedef uint128_t maxuint_t;

}

#elif defined(HAVE___INT128_T)

#define HAVE_INT128_T

#include <ostream>
#include <string>

namespace primecount {

typedef __int128_t int128_t;
typedef __int128_t maxint_t;

// uint128_t division is about 10% faster than int128_t division
// using GCC 4.8 on my Intel i7-4770 CPU from 2013
typedef __uint128_t uint128_t;
typedef __uint128_t maxuint_t;

inline std::ostream& operator<<(std::ostream& stream, uint128_t n)
{
  std::string str;

  while (n > 0)
  {
    str += '0' + n % 10;
    n /= 10;
  }

  if (str.empty())
    str = "0";

  stream << std::string(str.rbegin(), str.rend());
  return stream;
}

inline std::ostream& operator<<(std::ostream& stream, int128_t n)
{
  if (n < 0)
  {
    stream << "-";
    n = -n;
  }

  stream << (uint128_t) n;
  return stream;
}

} // namespace

#else /* int128_t not supported */

namespace primecount {

typedef int64_t maxint_t;
typedef int64_t maxuint_t;

}

#endif

#endif /* PTYPES_HPP */
