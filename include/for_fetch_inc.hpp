///
/// @file  for_fetch_inc.hpp
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef FOR_FETCH_INC
#define FOR_FETCH_INC

/// for (b = start; condition; b = atomic_i++)
/// This is a for loop with with dynamic thread scheduling for use
/// inside of a parallel region. We use this instead of OpenMP's
/// dynamic thread scheduling because the Clang compiler currently
/// (2021) has a severe scaling issue on PCs/servers with a large
/// number of CPU cores. GCC has no such scaling issue.
///
#define for_fetch_inc(atomic_i, start, condition) \
  _Pragma("omp single") \
  atomic_i = start; \
  for (auto b = atomic_i.fetch_add(1, memory_order_relaxed); \
       condition; b = atomic_i.fetch_add(1, memory_order_relaxed))

#endif
