///
/// @file  primecount-config.hpp
/// @brief Default CPU cache sizes and maximum CPU cache line size
///        that will be used by primecount's algorithms.
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PRIMECOUNT_CONFIG_HPP
#define PRIMECOUNT_CONFIG_HPP

/// Per CPU core L1 data cache size in bytes.
/// For big.LITTLE CPUs pick the smallest L1 data cache.
#ifndef L1D_CACHE_SIZE
  #define L1D_CACHE_SIZE (64 << 10)
#endif

/// Per CPU core L2 cache size in bytes.
/// For big.LITTLE CPUs pick the smallest L2 cache.
/// If your CPU's L2 cache is shared by multiple physical CPU cores
/// then L2_CACHE_SIZE should be set to (SHARED_L2_CACHE_SIZE /
/// number of physical CPU cores sharing the L2 cache).
#ifndef L2_CACHE_SIZE
  #define L2_CACHE_SIZE (512 << 10)
#endif

/// Maximum CPU cache line size in bytes (of all CPU types that
/// will be produced over the next few decades).
/// In order to prevent false sharing when using a mutex (or atomic
/// variable) we need to ensure that this mutex is stored on a
/// cache line where no other data is stored. We achieve this by
/// adding MAX_CACHE_LINE_SIZE bytes before and after the mutex.
#ifndef MAX_CACHE_LINE_SIZE
  #define MAX_CACHE_LINE_SIZE 512
#endif

#endif
