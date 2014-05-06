///
/// @file  Pk.hpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PK_HPP
#define PK_HPP

#include <primecount.hpp>
#include <stdint.h>

namespace primecount {

int64_t P2(int64_t x, int64_t a, int64_t y);
int64_t P3(int64_t x, int64_t a, int threads = MAX_THREADS);

} // namespace primecount

#endif
