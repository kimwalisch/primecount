///
/// @file  pi_deleglise_rivat.cpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>

#include <iostream>

namespace primecount {

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
///
int64_t pi_deleglise_rivat(int64_t x, int threads)
{
  return pi_deleglise_rivat_parallel1(x, threads);
}

} // namespace primecount
