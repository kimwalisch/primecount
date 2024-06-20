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
#include <stdint.h>

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
inline uint64_t get_xcr0()
{
#if defined(_MSC_VER)
  return _xgetbv(0);
#else
  uint32_t eax;
  uint32_t edx;

  __asm__ ("xgetbv" : "=a"(eax), "=d"(edx) : "c"(0));
  return eax | (uint64_t(edx) << 32);
#endif
}

inline bool run_cpuid_avx512_bmi2()
{
  int abcd[4];

  run_cpuid(1, 0, abcd);

  int osxsave_mask = (1 << 27);

  // Ensure OS supports extended processor state management
  if ((abcd[2] & osxsave_mask) != osxsave_mask)
    return false;

  uint64_t ymm_mask = XSTATE_SSE | XSTATE_YMM;
  uint64_t zmm_mask = XSTATE_SSE | XSTATE_YMM | XSTATE_ZMM;
  uint64_t xcr0 = get_xcr0();

  // Check AVX OS support
  if ((xcr0 & ymm_mask) != ymm_mask)
    return false;

  // Check AVX512 OS support
  if ((xcr0 & zmm_mask) != zmm_mask)
    return false;

  run_cpuid(7, 0, abcd);

  if ((abcd[1] & bit_BMI2) != bit_BMI2)
    return false;

  // AVX512F, AVX512VPOPCNTDQ
  return ((abcd[1] & bit_AVX512F) == bit_AVX512F &&
          (abcd[2] & bit_AVX512_VPOPCNTDQ) == bit_AVX512_VPOPCNTDQ);
}

/// Initialized at startup
bool cpu_supports_avx512_bmi2 = run_cpuid_avx512_bmi2();

} // namespace

#endif
