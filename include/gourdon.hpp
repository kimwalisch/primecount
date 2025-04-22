///
/// @file  gourdon.hpp
/// @brief Function declarations related to Xavier Gourdon's prime
///        counting function algorithm.
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef GOURDON_HPP
#define GOURDON_HPP

#include <cpu_arch_macros.hpp>
#include <int128_t.hpp>
#include <print.hpp>

#include <stdint.h>

namespace primecount {

int64_t pi_gourdon(int64_t x, int threads);
int64_t pi_gourdon_64(int64_t x, int threads, bool print = is_print());
int64_t Sigma(int64_t x, int64_t y, int threads, bool print = is_print());
int64_t Phi0(int64_t x, int64_t y, int64_t z, int64_t k, int threads, bool print = is_print());
int64_t AC(int64_t x, int64_t y, int64_t z, int64_t k, int threads, bool print = is_print());
int64_t B(int64_t x, int64_t y, int threads, bool print = is_print());
int64_t D(int64_t x, int64_t y, int64_t z, int64_t k, int64_t d_approx, int threads, bool print = is_print());

#if defined(ENABLE_PORTABLE_POPCNT64) || \
    defined(ENABLE_AVX512_VPOPCNT) || \
    defined(ENABLE_ARM_SVE)
  int64_t D_default(int64_t x, int64_t y, int64_t z, int64_t k, int64_t d_approx, int threads, bool print);
#endif

#if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  int64_t D_multiarch_avx512(int64_t x, int64_t y, int64_t z, int64_t k, int64_t d_approx, int threads, bool print);
#endif

#if defined(ENABLE_MULTIARCH_ARM_SVE)
  int64_t D_default(int64_t x, int64_t y, int64_t z, int64_t k, int64_t d_approx, int threads, bool print);
#endif

#ifdef HAVE_INT128_T

int128_t pi_gourdon(int128_t x, int threads);
int128_t pi_gourdon_128(int128_t x, int threads, bool print = is_print());
int128_t Sigma(int128_t x, int64_t y, int threads, bool print = is_print());
int128_t Phi0(int128_t x, int64_t y, int64_t z, int64_t k, int threads, bool print = is_print());
int128_t AC(int128_t x, int64_t y, int64_t z, int64_t k, int threads, bool print = is_print());
int128_t B(int128_t x, int64_t y, int threads, bool print = is_print());
int128_t D(int128_t x, int64_t y, int64_t z, int64_t k, int128_t d_approx, int threads, bool print = is_print());

#if defined(ENABLE_PORTABLE_POPCNT64) || \
    defined(ENABLE_AVX512_VPOPCNT) || \
    defined(ENABLE_ARM_SVE)
  int128_t D_default(int128_t x, int64_t y, int64_t z, int64_t k, int128_t d_approx, int threads, bool print);
#endif

#if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  int128_t D_multiarch_avx512(int128_t x, int64_t y, int64_t z, int64_t k, int128_t d_approx, int threads, bool print);
#endif

#if defined(ENABLE_MULTIARCH_ARM_SVE)
  int128_t D_multiarch_arm_sve(int128_t x, int64_t y, int64_t z, int64_t k, int128_t d_approx, int threads, bool print);
#endif

#endif

} // namespace

#endif
