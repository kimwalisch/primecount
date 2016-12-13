///
/// @file  S1.hpp
/// @brief S1 function declarations.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <int128_t.hpp>
#include <stdint.h>

namespace primecount {

int64_t S1(int64_t x,
           int64_t y,
           int64_t c,
           int threads);

#ifdef HAVE_INT128_T

int128_t S1(int128_t x,
            int64_t y,
            int64_t c,
            int threads);

#endif

} // namespace
