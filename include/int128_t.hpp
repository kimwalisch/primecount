///
/// @file   int128_t.hpp
/// @brief  Convenience functions and typedefs for the non standard
///         __int128_t integer type.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef INT128_T_HPP
#define INT128_T_HPP

#if defined(HAS_INT128_T)

#include <stdint.h>

#elif defined(HAS___INT128_T)

// int128_t is a typedef for __int128_t
#define HAS_INT128_T

#include <ostream>
#include <string>

namespace primecount {

typedef __int128_t int128_t;

std::ostream& operator<<(std::ostream& stream, int128_t n)
{
  std::string str;
  while (n > 0)
  {
    str += '0' + n % 10;
    n /= 10;
  }
  stream << std::string(str.rbegin(), str.rend());
  return stream;
}

} // namespace

#endif /* HAS___INT128_T */

#endif /* INT128_T_HPP */
