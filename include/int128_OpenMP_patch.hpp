///
/// @file   int128_OpenMP_patch.hpp
/// @brief  Add missing int128_t support for OpenMP multi-threading.
///         This is mainly required for LLVM/Clang on both Linux
///         and Windows. Without these pragmas there will often be
///         "unresolved external symbol" linker errors for symbols
///         like e.g. __atomic_load, __atomic_compare_exchange.
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
