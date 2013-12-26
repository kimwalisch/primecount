///
/// @file  PhiTiny.cpp
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "PhiTiny.h"
#include <stdint.h>

namespace primecount {

const int32_t PhiTiny::primes_[7] = { 0, 2, 3, 5, 7, 11, 13 };

/// prime_products_[n] = \prod_{i=1}^{n} primes_[i]
const int32_t PhiTiny::prime_products_[7] = { 1, 2, 6, 30, 210, 2310, 30030 };

/// totients_[n] = \prod_{i=1}^{n} (primes_[i] - 1)
const int32_t PhiTiny::totients_[7] = { 1, 1, 2, 8, 48, 480, 5760 };

} // namespace primecount
