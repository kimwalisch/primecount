///
/// @file  CPUID.hpp
/// @brief POPCNT detection fo x86 and x86-64 CPUs.
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef CPUID_HPP
#define CPUID_HPP

// Enable on x86 and x86-64 CPUs
#if defined(__x86_64__) || \
    defined(__i386__) || \
    defined(_M_X64) || \
    defined(_M_IX86)

// Check if CPUID POPCNT runtime check is needed
#if !(defined(__POPCNT__) || \
     (defined(_MSC_VER) && defined(__AVX__)))

#define HAS_CPUID_POPCNT

#if defined(_MSC_VER)
  #include <intrin.h>
#endif

namespace {

void run_CPUID(int eax, int ecx, int* abcd)
{
#if defined(_MSC_VER)
  __cpuidex(abcd, eax, ecx);
#else
  int ebx = 0;
  int edx = 0;

  #if defined(__i386__) && \
      defined(__PIC__)
    /* in case of PIC under 32-bit EBX cannot be clobbered */
    __asm__ ("movl %%ebx, %%edi;"
             "cpuid;"
             "xchgl %%ebx, %%edi;"
             : "=D" (ebx),
               "+a" (eax),
               "+c" (ecx),
               "=d" (edx));
  #else
    __asm__ ("cpuid;"
             : "+b" (ebx),
               "+a" (eax),
               "+c" (ecx),
               "=d" (edx));
  #endif

  abcd[0] = eax;
  abcd[1] = ebx;
  abcd[2] = ecx;
  abcd[3] = edx;
#endif
}

bool run_CPUID_POPCNT()
{
  // %ecx POPCNT bit flag
  int bit_POPCNT = 1 << 23;
  int abcd[4];

  run_CPUID(1, 0, abcd);
  return (abcd[2] & bit_POPCNT) == bit_POPCNT;
}

/// Initialized at startup
const bool CPUID_POPCNT = run_CPUID_POPCNT();

} // namespace

#endif // __POPCNT__
#endif // x86

#endif
