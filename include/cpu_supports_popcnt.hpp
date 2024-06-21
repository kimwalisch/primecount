///
/// @file  cpu_supports_popcnt.hpp
/// @brief POPCNT detection fo x86 and x86-64 CPUs.
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef CPU_SUPPORTS_POPCNT_HPP
#define CPU_SUPPORTS_POPCNT_HPP

// Enable CPUID for POPCNT on x86 and x86-64 CPUs.
// This is required because not all x86 and x86-64 CPUs
// support the POPCNT instruction.
#if defined(__x86_64__) || \
    defined(__i386__) || \
    defined(_M_X64) || \
    defined(_M_IX86)

// Both GCC and Clang (even Clang on Windows) define the __POPCNT__
// macro if the user compiles with -mpopcnt. The __POPCNT__
// macro is even defined if the user compiles with other flags
// such as -mavx or -march=native.
#if defined(__POPCNT__)
  #define HAS_POPCNT

// The MSVC compiler does not support a POPCNT macro, but if the user
// compiles with e.g. /arch:AVX or /arch:AVX512 then MSVC defines
// the __AVX__ macro and POPCNT is also supported.
#elif defined(_MSC_VER) && defined(__AVX__)
  #define HAS_POPCNT
#endif

#if !defined(HAS_POPCNT)
#define ENABLE_CPUID_POPCNT

namespace primecount {

bool has_cpuid_popcnt();

} // namespace

namespace {

/// Initialized at startup
bool cpu_supports_popcnt = primecount::has_cpuid_popcnt();

} // namespace

#endif // !defined(HAS_POPCNT)
#endif // x86 or x86-64

#endif
