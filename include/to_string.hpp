///
/// @file   to_string.hpp
/// @brief  to_string(n) implementation with 128-bit integer support.
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef TO_STRING_HPP
#define TO_STRING_HPP

#include <int128_t.hpp>
#include <string>

#if defined(HAVE_NON_STANDARD__INT128_T)

namespace primecount {

/// defined in util.cpp
std::string to_string(int128_t x);
std::string to_string(uint128_t x);

} // namespace

#else

namespace primecount {

using std::to_string;

} // namespace

#endif

#endif
