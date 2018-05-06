///
/// @file  pi_primesieve.cpp
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primesieve.hpp>

#include <stdint.h>

namespace primecount {

int64_t pi_primesieve(int64_t x)
{
  if (x < 2)
    return 0;
  else
    return primesieve::count_primes(0, x);
}

} // namespace
