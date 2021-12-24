///
/// @file   int128_t.hpp
/// @brief  Support for int128_t, uint128_t types.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef INT128_T_HPP
#define INT128_T_HPP

#include <stdint.h>
#include <string>

/// If INT128_MAX is defined we know that int128_t and
/// uint128_t are available in <stdint.h>.
///
#if defined(INT128_MAX) && \
   !defined(DISABLE_INT128)

#define HAVE_INT128_T

namespace primecount {

using maxint_t = int128_t;
using maxuint_t = uint128_t;
using std::to_string;

}

/// The __int128_t type (GCC/Clang) is not well supported by
/// the C++ standard library (in 2016) so we have to define
/// some functions ourselves. We also define typedefs so we
/// can use int128_t instead of __int128_t. Once this is done
/// int128_t can be used like a regular integer type.
///
#elif defined(__SIZEOF_INT128__) && \
     !defined(DISABLE_INT128)

#define HAVE_INT128_T
#define HAVE_NON_STANDARD__INT128_T

#include <ostream>

namespace primecount {

using int128_t = __int128_t;
using uint128_t = __uint128_t;
using maxint_t = __int128_t;
using maxuint_t = __uint128_t;

/// defined in util.cpp
std::string to_string(int128_t x);
std::string to_string(uint128_t x);
std::ostream& operator<<(std::ostream& stream, int128_t n);
std::ostream& operator<<(std::ostream& stream, uint128_t n);

} // namespace

#else // int128_t not supported

namespace primecount {

using maxint_t = int64_t;
using maxuint_t = uint64_t;
using std::to_string;

}

#endif

#endif
