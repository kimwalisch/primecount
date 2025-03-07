///
/// @file  cpu_arch_macros.hpp
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef CPU_ARCH_MACROS_HPP
#define CPU_ARCH_MACROS_HPP

// Needed for __has_include
#include <macros.hpp>

#if defined(__ARM_FEATURE_SVE) && \
    __has_include(<arm_sve.h>)
  #define ENABLE_ARM_SVE
#elif defined(__AVX512F__) && \
      defined(__AVX512VPOPCNTDQ__) && \
      __has_include(<immintrin.h>)
  #define ENABLE_AVX512_VPOPCNT
#elif defined(ENABLE_MULTIARCH_ARM_SVE)
  #define ENABLE_PORTABLE
#elif defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  #define ENABLE_PORTABLE
#else
  #define ENABLE_PORTABLE
#endif

#endif
