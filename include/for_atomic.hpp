///
/// @file  for_atomic.hpp
///        The for_atomic_inc() macro is a workaround for a severe
///        scaling issue in Clang's OpenMP library when using a parallel
///        for loop with schedule(dynamic). More details can be found
///        in my bug report: https://bugs.llvm.org/show_bug.cgi?id=49588
///        After having experienced the above scaling issue I don't
///        trust OpenMP dynamic scheduling anymore and I avoid using it.
///
///        Note that the atomic for loops in this file are optimal for
///        primecount's use case. Unlike GCC's OpenMP library which
///        uses atomics with memory_order_acq_rel our implementation
///        uses atomics with memory_order_relaxed because our threads
///        are completely independent from each other. This improves
///        performance on CPU architectures with a weak memory model
///        such as ARM64. Also since our algorithms perform best using
///        a chunk size of 1 we can hardcode that value.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef FOR_ATOMIC
#define FOR_ATOMIC

#include <atomic>

/// for_atomic_inc(start, condition)
/// Is the same as:
///
/// #pragma omp for nowait schedule(dynamic)
/// for (auto b = start; condition; b++)
///
#define for_atomic_inc(atomic_b, condition) \
  for (auto b = atomic_b.fetch_add(1, std::memory_order_relaxed); \
       condition; b = atomic_b.fetch_add(1, std::memory_order_relaxed))

/// parallel_for_atomic_inc(start, condition)
/// Is the same as:
///
/// #pragma omp parallel for schedule(dynamic) num_threads(threads) reduction(+: sum)
/// for (auto b = start; condition; b++)
///
#define parallel_for_atomic_inc(start, condition) \
  std::atomic<decltype(start)> atomic_b(start); \
  _Pragma("omp parallel num_threads(threads) reduction(+: sum)") \
  for (decltype(start) b = atomic_b.fetch_add(1, std::memory_order_relaxed); \
       condition; b = atomic_b.fetch_add(1, std::memory_order_relaxed))

#endif
