///
/// @file  pi.cpp
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <stdint.h>

namespace primecount {

/// Alias for fastest pi(x) implementation
int64_t pi(int64_t x, int threads /* = MAX_THREADS */)
{
  return pi_lehmer(x, threads);
}

} // namespace primecount
