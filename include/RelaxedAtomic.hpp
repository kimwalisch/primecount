///
/// @file  RelaxedAtomic.hpp
///        LLVM/OpenMP dynamic scheduling caused a severe scaling issue
///        in primecount: https://bugs.llvm.org/show_bug.cgi?id=49588
///        By default OpenMP dynamic scheduling may process iterations
///        in random order, it is likely that this causes many cache
///        misses when computing the easy special leaves in primecount.
///        It is possible to use schedule(monotonic:dynamic) to ensure
///        that iterations are processed in sequential order. This
///        seems to fix primecount's scaling issue.
///
///        However as a precaution measure we try to avoid OpenMP
///        dynamic scheduling in primecount and instead implement
///        schedule(monotonic:dynamic, 1) ourselves using relaxed
///        atomics.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef RELAXEDATOMIC_HPP
#define RELAXEDATOMIC_HPP

#include <atomic>

template <typename T>
class RelaxedAtomic
{
public:
  RelaxedAtomic(T n) : atomic_(n) { }
  /// Postfix Increment
  T operator++(int)
  {
    return atomic_.fetch_add(1, std::memory_order_relaxed);
  }
  char unused()
  {
    pad1[101] = 101;
    pad2[102] = 102;
    return pad1[101] | pad2[102];
  }
private:
  // Use padding to avoid CPU false sharing
  char pad1[512];
  std::atomic<T> atomic_;
  char pad2[512];
};

#endif
