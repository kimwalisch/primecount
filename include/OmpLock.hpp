///
/// @file   OmpLock.hpp
/// @brief  The OmpLock and LockGuard classes are RAII-style
///         wrappers for OpenMP locks.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef OMPLOCK_HPP
#define OMPLOCK_HPP

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
  OmpLock()
  {
    omp_init_lock(&lock_);
  }
  ~OmpLock()
  {
    omp_destroy_lock(&lock_);
  }
  char unused()
  {
    pad1[101] = 101;
    pad2[102] = 102;
    return pad1[101] | pad2[102];
  }

  // Use padding to avoid CPU false sharing
  char pad1[512];
  omp_lock_t lock_;
  char pad2[512];
};

class LockGuard
{
public:
  LockGuard(OmpLock& lock)
  {
    lock_ = &lock.lock_;
    omp_set_lock(lock_);
  }
  ~LockGuard()
  {
    omp_unset_lock(lock_);
  }

private:
  omp_lock_t* lock_;
};

} // namespace

#endif
