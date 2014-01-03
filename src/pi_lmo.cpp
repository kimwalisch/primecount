///
/// @file  pi_lmo.cpp
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <stdint.h>
#include <iostream>

namespace primecount {

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3)) operations, O(x^0.5) space.
/// @note O(x^0.5) space is due to parallel P2(x, a).
///
int64_t pi_lmo(int64_t x, int threads)
{
  std::cerr << "Error: pi_lmo(x) not yet implemented." << std::endl;
  return -1;
}

} // namespace primecount
