///
/// @file  cpuid.hpp
/// @brief CPUID for x86 and x86-64 CPUs.
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef CPUID_HPP
#define CPUID_HPP

#if defined(_MSC_VER)
  #include <intrin.h>
#endif

namespace {

inline void run_cpuid(int eax, int ecx, int* abcd)
{
#if defined(_MSC_VER)
  __cpuidex(abcd, eax, ecx);
#else
  int ebx = 0;
  int edx = 0;

  #if defined(__i386__) && \
      defined(__PIC__)
    // In case of PIC under 32-bit EBX cannot be clobbered
    __asm__ ("movl %%ebx, %%edi;"
             "cpuid;"
             "xchgl %%ebx, %%edi;"
             : "+a" (eax),
               "=D" (ebx),
               "+c" (ecx),
               "=d" (edx));
  #else
    __asm__ ("cpuid"
             : "+a" (eax),
               "+b" (ebx),
               "+c" (ecx),
               "=d" (edx));
  #endif

  abcd[0] = eax;
  abcd[1] = ebx;
  abcd[2] = ecx;
  abcd[3] = edx;
#endif
}

} // namespace

#endif
