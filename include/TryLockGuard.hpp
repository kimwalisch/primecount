 ///
/// @file   TryLockGuard.hpp
/// @brief  The TryLockGuard class is a RAII-style wrapper
///         for a fast lock using std::atomic.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef TRYLOCKGUARD_HPP
#define TRYLOCKGUARD_HPP

#include <atomic>

namespace {
 
/// For printing the status it is OK to use a non-blocking
/// userspace lock because printing the status is a non
/// essential operation and hence even if the OS preempts
/// the thread holding the lock it won't cause any deadlocks
/// or performance issues, it will only delay the status
/// output.
///
struct TryLockGuard
{
  public:
  TryLockGuard(std::atomic<bool>& lock)
      : lock_(&lock)
  {
      bool expected = false;
      is_locked_ = lock.compare_exchange_weak(expected, true, std::memory_order_acquire);
  }
  ~TryLockGuard()
  {
      if (is_locked_)
      lock_->store(false, std::memory_order_release);
  }
  bool owns_lock() const
  {
      return is_locked_;
  }
  private:
  std::atomic<bool>* lock_ = nullptr;
  bool is_locked_ = false;
};

} // namespace

#endif
