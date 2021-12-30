///
/// @file   int128_t.hpp
/// @brief  Support for int128_t, uint128_t types.
///         The code in this file uses only "old" C++ features from
///         C++03 and C++98 because it is used in primecount's main
///         CMakeLists.txt (CMake build script) and must work without
///         any special compiler flags with all C++ compilers.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef INT128_T_HPP
#define INT128_T_HPP

#include <stdint.h>

/// If INT128_MAX is defined we know that int128_t and
/// uint128_t are available in <stdint.h>.
///
#if defined(INT128_MAX) && \
   !defined(DISABLE_INT128)

#define HAVE_INT128_T

namespace primecount {

typedef int128_t maxint_t;
typedef uint128_t maxuint_t;

} // namespace

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

typedef __int128_t int128_t;
typedef __uint128_t uint128_t;
typedef __int128_t maxint_t;
typedef __uint128_t maxuint_t;

/// defined in util.cpp
std::ostream& operator<<(std::ostream& stream, int128_t n);
std::ostream& operator<<(std::ostream& stream, uint128_t n);

} // namespace

#else // int128_t not supported

namespace primecount {

typedef int64_t maxint_t;
typedef uint64_t maxuint_t;

} // namespace

#endif

#endif
