///
/// @file  pi_primesieve.cpp
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <primesieve.hpp>

#include <stdint.h>

namespace primecount {

int64_t pi_primesieve(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  primesieve::set_num_threads(threads);
  int64_t pix = primesieve::count_primes(0, x);
  
  return pix;
}

} // namespace
