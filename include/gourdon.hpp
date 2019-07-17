///
/// @file  gourdon.hpp
/// @brief Function declarations related to Xavier Gourdon's prime
///        counting function algorithm.
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <int128_t.hpp>
#include <stdint.h>

namespace primecount {

int64_t pi_gourdon1(int64_t x, int64_t y, int64_t z, int64_t k);

int64_t phi0(int64_t x,
             int64_t y,
             int64_t z,
             int64_t k,
             int threads);

#ifdef HAVE_INT128_T

int128_t phi0(int128_t x,
              int64_t y,
              int64_t z,
              int64_t k,
              int threads);

#endif

} // namespace
