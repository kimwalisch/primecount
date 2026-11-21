///
/// @file  RelaxedAtomic.hpp
///        Padded atomic counter with relaxed memory ordering. Threads
///        use it to process loop iterations in sequential order, like
///        OpenMP's schedule(monotonic:dynamic, 1). LLVM/OpenMP's
///        default dynamic scheduling, which may process iterations
///        in random order, caused a severe scaling issue in
///        primecount, likely due to many cache misses when computing
///        the easy special leaves:
///        https://bugs.llvm.org/show_bug.cgi?id=49588
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef RELAXEDATOMIC_HPP
#define RELAXEDATOMIC_HPP

#include <macros.hpp>
#include <primecount-config.hpp>

#include <atomic>

namespace {

template <typename T>
class RelaxedAtomic
{
public:
  RelaxedAtomic(T n) : atomic_(n) { }
  // Postfix Increment
  T operator++(int)
  {
    return atomic_.fetch_add(1, std::memory_order_relaxed);
  }
private:
  // Use padding to avoid CPU false sharing
  MAYBE_UNUSED char pad1[MAX_CACHE_LINE_SIZE];
  std::atomic<T> atomic_;
  MAYBE_UNUSED char pad2[MAX_CACHE_LINE_SIZE];
};

} // namespace

#endif
