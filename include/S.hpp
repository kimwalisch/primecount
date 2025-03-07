///
/// @file  S.hpp
/// @brief The S1 and S2 functions are part of the Deleglise-Rivat
///        prime counting algorithm.
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef S_HPP
#define S_HPP

#include <cpu_arch_macros.hpp>
#include <int128_t.hpp>
#include <print.hpp>

#include <stdint.h>

namespace primecount {

int64_t S1(int64_t x, int64_t y, int64_t c, int threads, bool print = is_print());
int64_t S2_trivial(int64_t x, int64_t y, int64_t z, int64_t c, int threads, bool print = is_print());
int64_t S2_easy(int64_t x, int64_t y, int64_t z, int64_t c, int threads, bool print = is_print());
int64_t S2_hard(int64_t x, int64_t y, int64_t z, int64_t c, int64_t s2_hard_approx, int threads, bool print = is_print());

#if defined(ENABLE_PORTABLE) || \
    defined(ENABLE_ARM_SVE) || \
    defined(ENABLE_AVX512_VPOPCNT)
  int64_t S2_hard_default(int64_t x, int64_t y, int64_t z, int64_t c, int64_t s2_hard_approx, int threads, bool print);
#endif

#if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  int64_t S2_hard_multiarch_avx512(int64_t x, int64_t y, int64_t z, int64_t c, int64_t s2_hard_approx, int threads, bool print);
#endif

#if defined(ENABLE_MULTIARCH_ARM_SVE)
  int64_t S2_hard_multiarch_arm_sve(int64_t x, int64_t y, int64_t z, int64_t c, int64_t s2_hard_approx, int threads, bool print);
#endif

#ifdef HAVE_INT128_T

int128_t S1(int128_t x, int64_t y, int64_t c, int threads, bool print = is_print());
int128_t S2_trivial(int128_t x, int64_t y, int64_t z, int64_t c, int threads, bool print = is_print());
int128_t S2_easy(int128_t x, int64_t y, int64_t z, int64_t c, int threads, bool print = is_print());
int128_t S2_hard(int128_t x, int64_t y, int64_t z, int64_t c, int128_t s2_hard_approx, int threads, bool print = is_print());

#if defined(ENABLE_PORTABLE) || \
    defined(ENABLE_ARM_SVE) || \
    defined(ENABLE_AVX512_VPOPCNT)
  int128_t S2_hard_default(int128_t x, int64_t y, int64_t z, int64_t c, int128_t s2_hard_approx, int threads, bool print);
#endif

#if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  int128_t S2_hard_multiarch_avx512(int128_t x, int64_t y, int64_t z, int64_t c, int128_t s2_hard_approx, int threads, bool print);
#endif

#if defined(ENABLE_MULTIARCH_ARM_SVE)
  int128_t S2_hard_multiarch_arm_sve(int128_t x, int64_t y, int64_t z, int64_t c, int128_t s2_hard_approx, int threads, bool print);
#endif

#endif

} // namespace

#endif
