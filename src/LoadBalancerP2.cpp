///
/// @file  LoadBalancerP2.cpp
/// @brief This load balancer assigns work to the threads in the
///        computation of the 2nd partial sieve function.
///        It is used by the P2(x, a) and B(x, y) functions.
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <LoadBalancerP2.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>
#include <min.hpp>

#include <stdint.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>

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
  // Using more chunks per thread improves load
  // balancing but also adds some overhead.
  int64_t chunks_per_thread = 8;
  min_thread_dist_ = 1 << 23;

  low_ = min(low_, sieve_limit_);
  int64_t dist = sieve_limit_ - low_;
  int threads1 = ideal_num_threads(dist, threads, min_thread_dist_);

  // This is a better approximation of the number of
  // threads used by the get_work() method for large
  // sieving distances. This formula applies if each
  // thread sieves a single segment of size
  // get_min_thread_dist(thread_low) and most of
  // these segments are > min_thread_dist_.
  double threads2 = (sieve_limit_ * 3.0) / get_min_thread_dist(sieve_limit_);

  threads_ = std::min(threads1, (int) threads2);
  threads_ = in_between(1, threads_, threads);
  thread_dist_ = dist / (threads_ * chunks_per_thread);
  thread_dist_ = max(min_thread_dist_, thread_dist_);

  lock_.init(threads_);
}

int LoadBalancerP2::get_threads() const
{
  return threads_;
}

/// Ensure that the thread initialization, i.e. the
/// computation of PrimePi(low), uses at most 10%
/// of the entire thread computation.
/// Since PrimePi(low) uses O(low^(2/3)/log(low)^2) time,
/// sieving a distance of n = low^(2/3) * 5 uses
/// O(n log log n) time, which is more than 10x more.
///
int64_t LoadBalancerP2::get_min_thread_dist(int64_t low) const
{
  double low13 = std::cbrt(low);
  int64_t n = (int64_t) (low13 * low13 * 5);
  return std::max(min_thread_dist_, n);
}

/// The thread needs to sieve [low, high[
bool LoadBalancerP2::get_work(int64_t& low, int64_t& high)
{
  LockGuard lockGuard(lock_);
  print_status();

  // Calculate the remaining sieving distance
  low_ = min(low_, sieve_limit_);
  int64_t dist = sieve_limit_ - low_;

  // When a single thread is used (and printing is
  // disabled) we can set thread_dist to the entire
  // sieving distance as load balancing is only
  // useful for multi-threading.
  if (threads_ == 1)
  {
    if (!is_print_)
      thread_dist_ = dist;
  }
  else
  {
    min_thread_dist_ = get_min_thread_dist(low_);
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

  return low < sieve_limit_;
}

void LoadBalancerP2::print_status()
{
  if (is_print_)
  {
    double time = get_time();
    double old = time_;
    double threshold = 0.1;

    if ((time - old) >= threshold)
    {
      time_ = time;
      std::cout << "\rStatus: " << std::fixed << std::setprecision(precision_)
                << get_percent(low_, sieve_limit_) << '%' << std::flush;
    }
  }
}

} // namespace
