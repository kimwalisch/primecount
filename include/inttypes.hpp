///
/// @file   inttypes.hpp
/// @brief  Integer types used in primecount: int128_t, uint128_t,
///         intfast64_t, intfast128_t, maxint_t, maxuint_t.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef INTTYPES_HPP
#define INTTYPES_HPP

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

#endif /* HAVE_INT128_T */

namespace primecount {

/// Fastest 64-bit integer type for division.
/// On most Intel CPUs before 2015 unsigned 64-bit division is about
/// 10 percent faster than signed division. It is likely that in a few
/// years signed and unsigned division will run equally fast.
///
typedef uint64_t intfast64_t;

#if defined(HAVE_INT128_T)

/// Fastest 128-bit integer type for division.
/// On the author's Intel Core-i7 4770 CPU from 2013 using uint128_t
/// instead of int128_t gives 10 percent better performance.
///
typedef uint128_t intfast128_t;

#endif

}

#endif /* INTTYPES_HPP */
