///
/// @file  pi_primesieve.cpp
/// @brief Count primes using the primesieve C/C++ library which
///        uses a highly optimized parallel implementation of the
///        segmented sieve of Eratosthenes.
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
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
