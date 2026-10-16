///
/// @file  cpu_supports_avx2.hpp
/// @brief Detect if the x86 CPU supports AVX2 and POPCNT.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef CPU_SUPPORTS_AVX2_HPP
#define CPU_SUPPORTS_AVX2_HPP

namespace primecount {

bool has_avx2();

} // namespace

namespace {

/// Initialized at startup
const bool cpu_supports_avx2 = primecount::has_avx2();

} // namespace

#endif
