///
/// @file  cpu_supports_arm_sve.hpp
///        Check if the CPU supports the ARM SVE instruction set.
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef CPU_SUPPORTS_ARM_SVE_HPP
#define CPU_SUPPORTS_ARM_SVE_HPP

#include <macros.hpp>

#if defined(__has_builtin) && \
    __has_builtin(__builtin_cpu_supports)

// GCC <= 8 had a bug where __builtin_cpu_supports("avx") would not
// check if AVX was supported by the operating system. Hence we
// only enable __builtin_cpu_supports() for GCC > 8.
// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=85100
#if defined(__GNUC__)
  #if __GNUC__ <= 8
    #define GCC_MISSING_OS_CHECK
  #endif
#endif

#if !defined(GCC_MISSING_OS_CHECK)

namespace {

/// Initialized at startup
bool cpu_supports_sve = __builtin_cpu_supports("sve");

} // namespace

#endif // !defined(GCC_MISSING_OS_CHECK)
#endif // __has_builtin(__builtin_cpu_supports)

#endif
