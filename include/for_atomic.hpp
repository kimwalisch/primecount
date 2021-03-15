///
/// @file  for_atomic.hpp
///        The for_atomic_inc() macro is a workaround for a severe
///        scaling issue in Clang's OpenMP library when using a parallel
///        for loop with schedule(dynamic). More details can be found
///        in my bug report: https://bugs.llvm.org/show_bug.cgi?id=49588
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

/// for_atomic_inc() is a for loop with dynamic thread scheduling for use
/// inside of an OpenMP parallel region. We use this instead of OpenMP's
/// dynamic thread scheduling because the Clang compiler currently
/// (Clang 11, 2021) has a severe scaling issue on PCs & servers with a
/// large number of CPU cores. The scaling issue occurred when computing
/// AC(x) with x >= 1e22.
///
/// for_atomic_inc(start, condition, atomic_b)
/// Is the same as:
///
/// #pragma omp for schedule(dynamic)
/// for (auto b = start; condition; b++)
///
#define for_atomic_inc(start, condition, atomic_b) \
  /* for_atomic_inc() is used in the computation of the AC formula */ \
  /* where the individual threads are completely independent of each */ \
  /* other, hence there is no communication between the threads and the */ \
  /* threads only read from memory but do not write to memory. */ \
  /* Because of these constraints we can use relaxed atomics, i.e. our */ \
  /* program will work correctly even if the CPU reorders memory reads */ \
  /* and writes before or after the atomic_b variable we use as loop */ \
  /* counter. Note that the fetch_add() instruction can never be */ \
  /* reordered before compare_exchange_strong() by the CPU because this */ \
  /* would change the program behavior in single-thread mode. */ \
  for (decltype(start) is_first_thread = -1, \
       b = (atomic_b.compare_exchange_strong(is_first_thread, start, std::memory_order_relaxed)) \
       ? atomic_b.fetch_add(1, std::memory_order_relaxed) \
       : atomic_b.fetch_add(1, std::memory_order_relaxed); \
       condition; b = atomic_b.fetch_add(1, std::memory_order_relaxed))

/// for_atomic_add() is a for loop with dynamic thread scheduling for use
/// inside of an OpenMP parallel region. We use this instead of OpenMP's
/// dynamic thread scheduling because the Clang compiler currently
/// (Clang 11, 2021) has a severe scaling issue on PCs & servers with a
/// large number of CPU cores. The scaling issue occurred when computing
/// AC(x) with x >= 1e22.
///
/// for_atomic_add(start, condition, atomic_b, inc)
/// Is the same as:
///
/// #pragma omp for schedule(dynamic)
/// for (auto b = start; condition; b = (atomic_b += inc))
///
#define for_atomic_add(start, condition, atomic_b, inc) \
  /* for_atomic_add() is used in the computation of the AC formula */ \
  /* where the individual threads are completely independent of each */ \
  /* other, hence there is no communication between the threads and the */ \
  /* threads only read from memory but do not write to memory. */ \
  /* Because of these constraints we can use relaxed atomics, i.e. our */ \
  /* program will work correctly even if the CPU reorders memory reads */ \
  /* and writes before or after the atomic_b variable we use as loop */ \
  /* counter. Note that the fetch_add() instruction can never be */ \
  /* reordered before compare_exchange_strong() by the CPU because this */ \
  /* would change the program behavior in single-thread mode. */ \
  for (decltype(start) is_first_thread = -1, \
       b = (atomic_b.compare_exchange_strong(is_first_thread, start, std::memory_order_relaxed)) \
       ? atomic_b.fetch_add(inc, std::memory_order_relaxed) \
       : atomic_b.fetch_add(inc, std::memory_order_relaxed); \
       condition; b = atomic_b.fetch_add(inc, std::memory_order_relaxed))

/// parallel_for_atomic_inc() is a parallel for loop with dynamic thread
/// scheduling. We use this instead of OpenMP's dynamic thread scheduling
/// because the Clang compiler currently (Clang 11, 2021) has a severe
/// scaling issue on PCs & servers with a large number of CPU cores. The
/// scaling issue occurred when computing S2_easy(x) with x >= 1e22.
///
/// parallel_for_atomic_inc(start, condition)
/// Is the same as:
///
/// #pragma omp parallel for schedule(dynamic) num_threads(threads) reduction(+: sum)
/// for (auto b = start; condition; b++)
///
#define parallel_for_atomic_inc(start, condition) \
  /* parallel_for_atomic_inc() is used in the computation of the S2_easy */ \
  /* formula where the individual threads are completely independent of */ \
  /* each  other, hence there is no communication between the threads and */ \
  /* the threads only read from memory but do not write to memory. */ \
  /* Because of these constraints we can use relaxed atomics, i.e. our */ \
  /* program will work correctly even if the CPU reorders memory reads */ \
  /* and writes before or after the atomic_b variable we use as loop */ \
  /* counter. Note that the fetch_add() instruction can never be */ \
  /* reordered before compare_exchange_strong() by the CPU because this */ \
  /* would change the program behavior in single-thread mode. */ \
  std::atomic<decltype(start)> atomic_b(start); \
  _Pragma("omp parallel num_threads(threads) reduction(+: sum)") \
  for (decltype(start) b = atomic_b.fetch_add(1, std::memory_order_relaxed); \
       condition; b = atomic_b.fetch_add(1, std::memory_order_relaxed))

#endif
