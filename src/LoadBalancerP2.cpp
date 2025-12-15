///
/// @file  LoadBalancerP2.cpp
/// @brief This load balancer assigns work to the threads in the
///        computation of the 2nd partial sieve function.
///        It is used by the P2(x, a) and B(x, y) functions.
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <LoadBalancerP2.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>
#include <min.hpp>
#include <print.hpp>

#include <stdint.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

namespace primecount {

/// We need to sieve [sqrt(x), sieve_limit[
LoadBalancerP2::LoadBalancerP2(maxint_t x,
                               int64_t sieve_limit,
                               int threads,
                               bool is_print) :
  low_(isqrt(x)),
  sieve_limit_(sieve_limit),
  precision_(get_status_precision(x)),
  is_print_(is_print)
{
  low_ = min(low_, sieve_limit_);
  int64_t dist = sieve_limit_ - low_;

  // These load balancing settings work well on my
  // dual-socket AMD EPYC 7642 server with 192 CPU cores.
  min_thread_dist_ = 1 << 23;
  int max_threads = (int) std::pow(sieve_limit_, 1 / 3.7);
  threads = std::min(threads, max_threads);
  threads_ = ideal_num_threads(dist, threads, min_thread_dist_);
  lock_.init(threads_);

  // Using more chunks per thread improves load
  // balancing but also adds some overhead.
  int64_t chunks_per_thread = 8;
  thread_dist_ = dist / (threads_ * chunks_per_thread);
  thread_dist_ = max(min_thread_dist_, thread_dist_);
}

int LoadBalancerP2::get_threads() const
{
  return threads_;
}

/// The thread needs to sieve [low, high[
bool LoadBalancerP2::get_work(int64_t& low, int64_t& high)
{
  std::string status;
  bool has_work;

  {
    LockGuard lockGuard(lock_);

    if (is_print_)
      status = get_status();

    // Calculate the remaining sieving distance
    low_ = min(low_, sieve_limit_);
    int64_t dist = sieve_limit_ - low_;

    // When a single thread is used (and printing is disabled)
    // we can set thread_dist to the entire sieving distance
    // as load balancing is only useful for multi-threading.
    if (threads_ == 1)
    {
      if (!is_print_)
        thread_dist_ = dist;
    }
    else
    {
      // Ensure that the thread initialization i.e. the calculation
      // of PrimePi(low) uses less time than the actual computation.
      // Computing PrimePi(low) uses O(low^(2/3) / log(low)^2) time but
      // sieving a distance of n = low^(2/3) uses O(n log log n) time,
      // hence the sieving time is larger than the initialization time.
      // Using these settings for low = 2e14 the sieving time is 5x
      // larger than the initialization time on my AMD EPYC 7642 CPU.
      double low13 = std::cbrt(low_);
      int64_t low23 = (int64_t) (low13 * low13);
      min_thread_dist_ = std::max(min_thread_dist_, low23);
      thread_dist_ = max(min_thread_dist_, thread_dist_);

      // Reduce the thread distance near to end to keep all
      // threads busy until the computation finishes.
      int64_t max_thread_dist = dist / threads_;
      if (thread_dist_ > max_thread_dist)
        thread_dist_ = max(min_thread_dist_, max_thread_dist);
    }

    low = low_;
    low_ += thread_dist_;
    low_ = min(low_, sieve_limit_);
    high = low_;
    has_work = low < sieve_limit_;
  }

  // Printing to the terminal incurs a system call
  // and may hence be slow. Therefore, we do it
  // after having released the mutex.
  if (!status.empty())
    std::cout << status << std::flush;

  return has_work;
}

std::string LoadBalancerP2::get_status()
{
  double time = get_time();
  double old = time_;
  double threshold = 0.1;

  if ((time - old) >= threshold)
  {
    time_ = time;
    double percent = get_percent(low_, sieve_limit_);
    std::string status;
    status.reserve(40);

    // Clear the previous status line since multiple
    // threads may print the status out of order.
    // Max status string size: "Status: 100.12345%",
    // hence we clear using 18 space characters.

    status = "\r                  \rStatus: ";
    status += to_string(percent, precision_);
    status += '%';
    return status;
  }

  return std::string();
}

} // namespace
