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

inline void run_cpuid(int* eax, int* ebx, int* ecx, int* edx)
{
#if defined(_MSC_VER)

  int abcd[4];
  __cpuidex(abcd, *eax, *ecx);

  *eax = abcd[0];
  *ebx = abcd[1];
  *ecx = abcd[2];
  *edx = abcd[3];

#elif defined(__i386__) && \
      defined(__PIC__)

  // in case of PIC under 32-bit EBX cannot be clobbered
  __asm__ (
    "movl %%ebx, %%edi;"
    "cpuid;"
    "xchgl %%ebx, %%edi;"
    : "=D" (*ebx),
      "+a" (*eax),
      "+c" (*ecx),
      "=d" (*edx)
  );

#else

  __asm__ (
    "cpuid"
    : "=a" (*eax),
      "=b" (*ebx),
      "=c" (*ecx),
      "=d" (*edx)
    : "a" (*eax), "c" (*ecx)
  );

#endif
}

} // namespace

#endif
