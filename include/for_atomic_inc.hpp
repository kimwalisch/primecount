///
/// @file  for_atomic_inc.hpp
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef FOR_ATOMIC_INC
#define FOR_ATOMIC_INC

#include <atomic>

/// for_atomic_inc() is a for loop with dynamic thread scheduling for
/// use inside of an OpenMP parallel region. We use this instead of
/// OpenMP's dynamic thread scheduling because the Clang compiler
/// currently (Clang 11, 2021) has a severe scaling issue on PCs &
/// servers with a large number of CPU cores. The scaling issue
/// occurred when computing AC(x) with x >= 1e22. GCC has no such
/// scaling issue but for_atomic_inc() runs slightly faster on GCC as
/// well.
///
/// for_atomic_inc(start, condition, atomic_i)
/// Is the same as:
///
/// #pragma omp for schedule(dynamic)
/// for (auto b = start; condition; b++)
///
#define for_atomic_inc(start, condition, atomic_i) \
  _Pragma("omp single") \
  atomic_i = start; \
  for (auto b = atomic_i.fetch_add(1, std::memory_order_relaxed); \
       condition; b = atomic_i.fetch_add(1, std::memory_order_relaxed))

#endif
