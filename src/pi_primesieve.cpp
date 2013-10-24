///
/// @file  pi_primesieve.cpp
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.h>
#include <primesieve.hpp>
#include <stdint.h>

namespace primecount {

int64_t pi_primesieve(int64_t x, int threads /* = MAX_THREADS */)
{
  if (x < 2)
    return 0;

  return primesieve::parallel_count_primes(0, x, threads);
}

} // namespace primecount
