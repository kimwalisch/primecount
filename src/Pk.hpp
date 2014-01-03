///
/// @file  Pk.hpp
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PK_PRIMECOUNT_H
#define PK_PRIMECOUNT_H

#include <primecount.hpp>
#include <stdint.h>

namespace primecount {

int64_t P2(int64_t x, int64_t a, int threads = MAX_THREADS);
int64_t P3(int64_t x, int64_t a, int threads = MAX_THREADS);

} // namespace primecount

#endif
