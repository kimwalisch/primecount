///
/// @file  cpu_arch_macros.hpp
/// @brief Enable the fastest SIMD instruction set supported by the
///        compiler flags. Runtime dispatching (ENABLE_MULTIARCH_*)
///        is only kept for instruction sets that are even faster.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef CPU_ARCH_MACROS_HPP
#define CPU_ARCH_MACROS_HPP

#include <macros.hpp>

#if defined(__AVX512F__) && \
    defined(__AVX512BW__) && \
    defined(__AVX512VL__) && \
    defined(__AVX512VPOPCNTDQ__) && \
    defined(__POPCNT__) && \
    !defined(__i386__) && \
    __has_include(<immintrin.h>)
  #define ENABLE_AVX512_VPOPCNT
  #undef ENABLE_MULTIARCH_AVX512_VPOPCNT
  #undef ENABLE_MULTIARCH_AVX2

#elif defined(__AVX2__) && \
      defined(__POPCNT__) && \
      !defined(__i386__) && \
      __has_include(<immintrin.h>)
  #define ENABLE_AVX2
  #undef ENABLE_MULTIARCH_AVX2

// MSVC does not define __POPCNT__, but
// /arch:AVX2 also enables POPCNT.
#elif defined(__AVX2__) && \
      defined(_MSC_VER) && \
      !defined(_M_IX86) && \
      __has_include(<immintrin.h>)
  #define ENABLE_AVX2
  #undef ENABLE_MULTIARCH_AVX2

#elif defined(__ARM_FEATURE_SVE) && \
      __has_include(<arm_sve.h>)
  #define ENABLE_ARM_SVE
  #undef ENABLE_MULTIARCH_ARM_SVE

// ARM NEON is supported by all 64-bit ARM CPUs. 32-bit ARM
// uses the portable default algorithms, as NEON is optional
// on 32-bit ARM and we use AArch64-only NEON intrinsics.
#elif (defined(__aarch64__) || \
       defined(_M_ARM64)) && \
      __has_include(<arm_neon.h>)
  #define ENABLE_ARM_NEON
#endif

#endif
