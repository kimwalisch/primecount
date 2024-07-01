///
/// @file  phi_vector.hpp
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PHI_VECTOR_HPP
#define PHI_VECTOR_HPP

#include <PiTable.hpp>
#include <Vector.hpp>
#include <stdint.h>

namespace primecount {

/// Returns a vector with phi(x, i - 1) values such that
/// phi[i] = phi(x, i - 1) for 1 <= i <= a.
/// phi(x, a) counts the numbers <= x that are not
/// divisible by any of the first a primes.
///
Vector<int64_t> phi_vector(int64_t x,
                           int64_t a,
                           const Vector<uint32_t>& primes,
                           const PiTable& pi);

/// Returns a vector with phi(x, i - 1) values such that
/// phi[i] = phi(x, i - 1) for 1 <= i <= a.
/// phi(x, a) counts the numbers <= x that are not
/// divisible by any of the first a primes.
///
Vector<int64_t> phi_vector(int64_t x,
                           int64_t a,
                           const Vector<int64_t>& primes,
                           const PiTable& pi);

} // namespace

#endif
