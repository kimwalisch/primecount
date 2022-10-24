///
/// @file   OmpLock.hpp
/// @brief  The OmpLock and LockGuard classes are RAII-style
///         wrappers for OpenMP locks.
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef OMPLOCK_HPP
#define OMPLOCK_HPP

#include <macros.hpp>

#ifndef MAX_CACHE_LINE_SIZE
  #define MAX_CACHE_LINE_SIZE 512
#endif

#if defined(_OPENMP)
  #include <omp.h>
#else

// If OpenMP is disabled we define the functions used by
// the OmpLock and LockGuard classes as no-op.
namespace {

using omp_lock_t = int;

inline void omp_init_lock(omp_lock_t*) { }
inline void omp_destroy_lock(omp_lock_t*) { }
inline void omp_set_lock(omp_lock_t*) { }
inline void omp_unset_lock(omp_lock_t*) { }

} // namespace

#endif

namespace primecount {

struct OmpLock
{
  void init(int threads)
  {
    ASSERT(threads > 0);
    threads_ = threads;

    if (threads_ > 1)
      omp_init_lock(&lock_);
  }

  ~OmpLock()
  {
    if (threads_ > 1)
      omp_destroy_lock(&lock_);
  }

  bool is_initialized() const
  {
    return threads_ > 0;
  }

  // 0 = uninitialized lock
  int threads_ = 0;

  // Use padding to avoid CPU false sharing
  MAYBE_UNUSED char pad1[MAX_CACHE_LINE_SIZE];
  omp_lock_t lock_;
  MAYBE_UNUSED char pad2[MAX_CACHE_LINE_SIZE];
};

class LockGuard
{
public:
  LockGuard(OmpLock& lock)
  {
    ASSERT(lock.is_initialized());

    if (lock.threads_ > 1)
    {
      lock_ = &lock.lock_;
      omp_set_lock(lock_);
    }
  }

  ~LockGuard()
  {
    if (lock_)
      omp_unset_lock(lock_);
  }

private:
  omp_lock_t* lock_ = nullptr;
};

} // namespace

#endif
