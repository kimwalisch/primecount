///
/// @file   int128_OpenMP_patch.hpp
/// @brief  Add missing int128_t support for OpenMP multi-threading.
///         This is required if libatomic does not support 128-bit
///         integers e.g. for LLVM/Clang on Windows. When including
///         this header OpenMP will use critical sections instead of
///         atomics for summing 128-bit integers. Critical sections
///         are less performant than atomics, hence this should only
///         be used as a last resort if nothing else works.
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef INT128_OPENMP_PATCH_HPP
#define INT128_OPENMP_PATCH_HPP

// Requires OpenMP >= 4.0
#pragma omp declare reduction(+: primecount::int128_t : omp_out += omp_in)

#pragma omp declare reduction(+: primecount::uint128_t : omp_out += omp_in)

#endif
