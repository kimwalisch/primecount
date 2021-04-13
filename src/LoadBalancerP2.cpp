///
/// @file  LoadBalancerP2.cpp
/// @brief This load balancer assigns work to the threads in the
///        computation of the 2nd partial sieve function.
///        It is used by the P2(x, a) and B(x, y) functions.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <LoadBalancerP2.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <iomanip>

using namespace std;

namespace primecount {

LoadBalancerP2::LoadBalancerP2(maxint_t x,
                               int64_t sieve_limit,
                               int threads,
                               bool is_print) :
  low_(isqrt(x)),
  sieve_limit_(sieve_limit),
  precision_(get_status_precision(x)),
  is_print_(is_print)
{
  int64_t min_dist = 1 << 22;
  int64_t dist = sieve_limit_ - min(low_, sieve_limit_);
  thread_dist_ = dist / (threads * 8);
  thread_dist_ = max(min_dist, thread_dist_);
  threads_ = ideal_num_threads(threads, dist, thread_dist_);

  if (!is_print && threads_ == 1)
    thread_dist_ = dist;
}

int LoadBalancerP2::get_threads() const
{
  return threads_;
}

bool LoadBalancerP2::get_work(int64_t& low, int64_t& high)
{
  LockGuard lockGuard(lock_);

  print_status();
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
    cout << "\rStatus: " << fixed << setprecision(precision_)
         << get_percent(low_, sieve_limit_) << '%' << flush;
  }
}

} // namespace
