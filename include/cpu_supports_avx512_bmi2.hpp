///
/// @file  cpu_supports_avx512_bmi2.hpp
/// @brief Detect if the x86 CPU supports AVX512 and BMI2.
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef CPU_SUPPORTS_AVX512_BMI2_HPP
#define CPU_SUPPORTS_AVX512_BMI2_HPP

#include <cpuid.hpp>

#if defined(_MSC_VER)
  #include <immintrin.h>
#endif

// CPUID bits documentation:
// https://en.wikipedia.org/wiki/CPUID

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
  __asm__ __volatile__ ("xgetbv" : "=a" (xcr0) : "c" (0) : "%edx" );
#endif

  return xcr0;
}

inline bool run_cpuid_avx512_bmi2()
{
  int eax = 1;
  int ebx = 0;
  int ecx = 0;
  int edx = 0;

  run_cpuid(&eax, &ebx, &ecx, &edx);

  int osxsave_mask = (1 << 27);

  // Ensure OS supports extended processor state management
  if ((ecx & osxsave_mask) != osxsave_mask)
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

  eax = 7;
  ebx = 0;
  ecx = 0;
  edx = 0;

  run_cpuid(&eax, &ebx, &ecx, &edx);

  if ((ebx & bit_BMI2) != bit_BMI2)
    return false;

  // AVX512F, AVX512VPOPCNTDQ
  return ((ebx & bit_AVX512F) == bit_AVX512F &&
          (ecx & bit_AVX512_VPOPCNTDQ) == bit_AVX512_VPOPCNTDQ);
}

/// Initialized at startup
bool cpu_supports_avx512_bmi2 = run_cpuid_avx512_bmi2();

} // namespace

#endif
