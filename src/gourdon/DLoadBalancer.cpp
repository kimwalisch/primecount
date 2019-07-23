///
/// @file  DLoadBalancer.cpp
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "DLoadBalancer.hpp"
#include <primecount-internal.hpp>
#include <Sieve.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <print.hpp>

#include <stdint.h>
#include <cmath>

using namespace std;

namespace primecount {

DLoadBalancer::DLoadBalancer(maxint_t x, 
                             int64_t y,
                             int64_t z) :
  sum_(0),
  low_(0),
  max_low_(0),
  xz_((int64_t)(x / z)),
  segments_(1),
  time_(get_time())
{
  init_size();
  maxint_t x16 = iroot<6>(x);
  double alpha_y = get_alpha_y(x, y);
  smallest_hard_leaf_ = (int64_t) (x / (y * sqrt(alpha_y) * x16));
}

void DLoadBalancer::init_size()
{
  // start with a tiny segment_size as most
  // special leaves are in the first few segments
  // and we need to ensure that all threads are
  // assigned an equal amount of work
  int64_t sqrt_xz = isqrt(xz_);
  int64_t log = ilog(sqrt_xz);
  log = max(log, 1);
  segment_size_ = sqrt_xz / log;

  int64_t min_size = 1 << 9;
  segment_size_ = max(segment_size_, min_size);
  segment_size_ = Sieve::get_segment_size(segment_size_);

  // try to use a segment size that fits exactly
  // into the CPUs L1 data cache
  int64_t l1_dcache_size = 1 << 15;
  max_size_ = l1_dcache_size * 30;
  max_size_ = max(max_size_, sqrt_xz);
  max_size_ = Sieve::get_segment_size(max_size_);
}

maxint_t DLoadBalancer::get_sum() const
{
  return sum_;
}

bool DLoadBalancer::get_work(int64_t* low,
                             int64_t* segments,
                             int64_t* segment_size,
                             maxint_t sum,
                             Runtime& runtime)
{
  #pragma omp critical (get_work)
  {
    sum_ += sum;

    update(low, segments, runtime);

    *low = low_;
    *segments = segments_;
    *segment_size = segment_size_;
    low_ += segments_ * segment_size_;
  }

  return *low <= xz_;
}

void DLoadBalancer::update(int64_t* low,
                           int64_t* segments,
                           Runtime& runtime)
{
  if (*low > max_low_)
  {
    max_low_ = *low;
    segments_ = *segments;

    if (segment_size_ < max_size_)
      segment_size_ = min(segment_size_ * 2, max_size_);
    else
      update_segments(runtime);
  }

  // Most hard special leaves are located just past
  // smallest_hard_leaf_. In order to prevent assigning
  // the bulk of work to a single thread we reduce
  // the number of segments to a minimum.

  int64_t high = low_ + segments_ * segment_size_;

  if (smallest_hard_leaf_ >= low_ &&
      smallest_hard_leaf_ <= high)
  {
    segments_ = 1;
  }
}

/// Increase or decrease the number of segments based on the
/// remaining runtime. Near the end it is important that
/// threads run only for a short amount of time in order to
/// ensure all threads finish nearly at the same time.
///
void DLoadBalancer::update_segments(Runtime& runtime)
{
  double min_secs = 0.01;

  // Calculate remaining time
  double quot = (double) low_ / xz_;
  quot = std::max(min_secs, quot);
  double rem_secs = (get_time() - time_) * (1 / quot);
  double threshold = rem_secs / 4;

  // Each thread should run at least 10x
  // longer than its initialization time
  threshold = max(threshold, runtime.init * 10);
  threshold = max(threshold, min_secs);

  // divider must not be 0
  double divider = max(runtime.secs, min_secs / 10);
  double factor = threshold / divider;

  // Reduce the thread runtime if it is much
  // larger than its initialization time
  if (runtime.secs > min_secs &&
      runtime.secs > runtime.init * 1000)
  {
    double old = factor;
    factor = (runtime.init * 1000) / runtime.secs;
    factor = min(factor, old);
  }

  factor = in_between(0.5, factor, 2.0);
  double new_segments = round(segments_ * factor);
  segments_ = (int64_t) new_segments;
  segments_ = max(segments_, 1);
}

} // namespace
