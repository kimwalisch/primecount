///
/// @file  pi_lmo.cpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "internal.hpp"

#include <primecount.hpp>
#include <iostream>

namespace primecount {

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
///
int64_t pi_lmo(int64_t x, int threads)
{
  return pi_lmo_parallel1(x, threads);
}

} // namespace primecount
