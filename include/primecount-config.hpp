///
/// @file  primecount-config.hpp
/// @brief Default CPU cache sizes and maximum CPU cache line size
///        that will be used by primecount's algorithms.
///
///        There are CPUs with small caches and CPUs with large
///        caches. In order to compile a primecount binary that will
///        perform well on a wide variety of CPUs, it is recommended
///        to set L1_CACHE_SIZE and L2_CACHE_SIZE to values found in
///        CPUs with small or medium sized caches.
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PRIMECOUNT_CONFIG_HPP
#define PRIMECOUNT_CONFIG_HPP

#ifndef L1_CACHE_SIZE
  /// L1 data cache size in bytes (per CPU core).
  /// For big.LITTLE CPUs pick the smallest L1 data cache.
  #define L1_CACHE_SIZE (64 << 10)
#endif

#ifndef L2_CACHE_SIZE
  /// L2 cache size in bytes (per CPU core).
  /// For big.LITTLE CPUs pick the smallest L2 cache.
  /// If your CPU's L2 cache is shared by multiple physical CPU cores
  /// then L2_CACHE_SIZE should be set to (SHARED_L2_CACHE_SIZE /
  /// number of physical CPU cores sharing the L2 cache).
  #define L2_CACHE_SIZE (512 << 10)
#endif

#ifndef MAX_CACHE_LINE_SIZE
  /// Maximum CPU cache line size in bytes (of all CPU types that
  /// will be produced over the next few decades).
  /// In order to prevent false sharing when using a mutex (or atomic
  /// variable) we need to ensure that this mutex is stored on a
  /// cache line where no other data is stored. We achieve this by
  /// adding MAX_CACHE_LINE_SIZE bytes before and after the mutex.
  #define MAX_CACHE_LINE_SIZE 512
#endif

#endif
