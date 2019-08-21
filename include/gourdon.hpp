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

int64_t pi_gourdon(int64_t x, int threads);
int64_t pi_gourdon_64(int64_t x, int threads);

int64_t Sigma(int64_t x, int64_t y, int threads);
int64_t Phi0(int64_t x, int64_t y, int64_t z, int64_t k, int threads);
int64_t AC(int64_t x, int64_t y, int64_t z, int64_t k, int threads);
int64_t B(int64_t x, int64_t y, int threads);
int64_t D(int64_t x, int64_t y, int64_t z, int64_t k, int64_t d_approx, int threads);

#ifdef HAVE_INT128_T

int128_t pi_gourdon(int128_t x, int threads);
int128_t pi_gourdon_128(int128_t x, int threads);

int128_t Sigma(int128_t x, int64_t y, int threads);
int128_t Phi0(int128_t x, int64_t y, int64_t z, int64_t k, int threads);
int128_t AC(int128_t x, int64_t y, int64_t z, int64_t k, int threads);
int128_t B(int128_t x, int64_t y, int threads);
int128_t D(int128_t x, int64_t y, int64_t z, int64_t k, int128_t d_approx, int threads);

#endif

} // namespace
