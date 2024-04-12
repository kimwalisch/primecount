///
/// @file  CPUID_AVX512_BMI2.hpp
/// @brief Detect if the x86 CPU supports AVX512 and BMI2.
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef CPUID_AVX512_BMI2_HPP
#define CPUID_AVX512_BMI2_HPP

#include <CPUID.hpp>

#if defined(_MSC_VER)
  #include <intrin.h>
#endif

// %ebx bit flags
#define bit_BMI2    (1 << 8)
#define bit_AVX512F (1 << 16)

// %ecx bit flags
#define bit_AVX512_VPOPCNTDQ (1 << 14)

// xgetbv bit flags
#define XSTATE_SSE (1 << 1)
#define XSTATE_YMM (1 << 2)
#define XSTATE_ZMM (7 << 5)

namespace {

// Get Value of Extended Control Register
inline int get_xcr0()
{
  int xcr0;

#if defined(_MSC_VER)
  xcr0 = (int) _xgetbv(0);
#else
  __asm__ ("xgetbv" : "=a" (xcr0) : "c" (0) : "%edx" );
#endif

  return xcr0;
}

inline bool get_AVX512_BMI2()
{
  int abcd[4];

  run_CPUID(1, 0, abcd);

  int osxsave_mask = (1 << 27);

  // Ensure OS supports extended processor state management
  if ((abcd[2] & osxsave_mask) != osxsave_mask)
    return false;

  int ymm_mask = XSTATE_SSE | XSTATE_YMM;
  int zmm_mask = XSTATE_SSE | XSTATE_YMM | XSTATE_ZMM;

  int xcr0 = get_xcr0();

  // Check AVX OS support
  if ((xcr0 & ymm_mask) != ymm_mask)
    return false;

  // Check AVX512 OS support
  if ((xcr0 & zmm_mask) != zmm_mask)
    return false;

  run_CPUID(7, 0, abcd);

  if ((abcd[1] & bit_BMI2) != bit_BMI2)
    return false;

  // AVX512F, AVX512VPOPCNTDQ
  return ((abcd[1] & bit_AVX512F) == bit_AVX512F &&
          (abcd[2] & bit_AVX512_VPOPCNTDQ) == bit_AVX512_VPOPCNTDQ);
}

/// Initialized at startup
const bool cpuid_AVX512_BMI2 = get_CPUID_AVX512_BMI2();

} // namespace

#endif
