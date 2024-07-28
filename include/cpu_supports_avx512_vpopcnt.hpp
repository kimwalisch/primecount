///
/// @file  cpu_supports_avx512_vpopcnt.hpp
/// @brief Detect if the x86 CPU supports AVX512 and BMI2.
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef CPU_SUPPORTS_AVX512_VPOPCNT_HPP
#define CPU_SUPPORTS_AVX512_VPOPCNT_HPP

namespace primecount {

bool has_cpuid_avx512_vpopcnt();

} // namespace

namespace {

/// Initialized at startup
const bool cpu_supports_avx512_vpopcnt = primecount::has_cpuid_avx512_vpopcnt();

} // namespace

#endif
