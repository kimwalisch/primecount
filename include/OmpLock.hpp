///
/// @file   OmpLock.hpp
/// @brief  The OmpLock and TryLock classes are RAII-style
///         wrappers for OpenMP locks.
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef OMPLOCK_HPP
#define OMPLOCK_HPP

#include <omp.h>

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
  omp_lock_t lock_;
};

class TryLock
{
public:
  TryLock(OmpLock& lock)
  {
    lock_ = &lock.lock_;
    ownsLock_ = (omp_test_lock(lock_) != 0);
  }
  ~TryLock()
  {
    if (ownsLock_)
      omp_unset_lock(lock_);
  }
  bool ownsLock() const
  {
    return ownsLock_;
  }

private:
  omp_lock_t* lock_;
  bool ownsLock_;
};

} // namespace

#endif
