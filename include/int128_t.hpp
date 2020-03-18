///
/// @file   int128_t.hpp
/// @brief  Support for int128_t, uint128_t types.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef INT128_T_HPP
#define INT128_T_HPP

#include <stdint.h>

/// If INT128_MAX is defined we know that int128_t and
/// uint128_t are available in stdint.h.
///
#if defined(INT128_MAX)
#define HAVE_INT128_T

namespace primecount {

using maxint_t = int128_t;
using maxuint_t = uint128_t;

}

/// The __int128_t type (GCC/Clang) is not well supported by
/// the C++ standard library (in 2016) so we have to define
/// some functions ourselves. We also define typedefs so we
/// can use int128_t instead of __int128_t. Once this is done
/// int128_t can be used like a regular integer type.
///
#elif defined(__SIZEOF_INT128__)
#define HAVE_INT128_T

#include <ostream>
#include <string>

namespace primecount {

using int128_t = __int128_t;
using uint128_t = __uint128_t;
using maxint_t = __int128_t;
using maxuint_t = __uint128_t;

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

#else // int128_t not supported

namespace primecount {

using maxint_t = int64_t;
using maxuint_t = uint64_t;

}

#endif

#endif
