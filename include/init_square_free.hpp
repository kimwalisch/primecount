///
/// @file  init_square_free.hpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef INIT_SQUARE_FREE_HPP
#define INIT_SQUARE_FREE_HPP

#include <stdint.h>
#include <vector>

namespace primecount {

/// Generate vectors containing n values which satisfy:
/// is_square_free(n) && && !is_prime(n) && primes[i] < least_prime_factor[n].
///
void init_square_free_candidates(std::vector<std::vector<int32_t> >*,
                                 std::vector<int32_t>&,
                                 std::vector<int32_t>&,
                                 std::vector<int32_t>&,
                                 std::vector<int32_t>&,
                                 int64_t c,
                                 int64_t y);

/// Initialize the square free iterators.
/// This version is for use in a single-threaded implementation.
///
void init_square_free_iters(std::vector<std::vector<int32_t>::iterator >*,
                            std::vector<std::vector<int32_t> >&);

/// Initialize the square free iterators.
/// This version is for use in a parallel implementation.
///
void init_square_free_iters(std::vector<std::vector<int32_t>::iterator >*,
                            std::vector<std::vector<int32_t> >&,
                            std::vector<int32_t>&,
                            int64_t c,
                            int64_t x,
                            int64_t y,
                            int64_t low);

} // namespace

#endif
