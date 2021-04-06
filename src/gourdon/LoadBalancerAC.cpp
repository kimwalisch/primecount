///
/// @file  LoadBalancerAC.cpp
/// @brief This load balancer assigns work to the threads in the
///        computation of the A & C formulas (AC.cpp) in
///        Xavier Gourdon's algorithm.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <LoadBalancerAC.hpp>
#include <SegmentedPiTable.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <algorithm>

namespace {

// CPU cache sizes per core
const int64_t l1_cache_size = 32 << 10;
const int64_t l2_cache_size = 256 << 10;

using primecount::SegmentedPiTable;
const int64_t numbers_per_byte = SegmentedPiTable::numbers_per_byte();

// Minimum segment size = 1 KiB
const int64_t min_segment_size = (1 << 10) * numbers_per_byte;

} // namespace

namespace primecount {

LoadBalancerAC::LoadBalancerAC(int64_t sqrtx,
                               int64_t y,
                               int threads) :
  low_(0),
  sqrtx_(sqrtx),
  x14_(isqrt(sqrtx)),
  y_(y),
  threads_(threads)
{
  if (threads_ == 1)
    segment_size_ = std::max(x14_, l2_cache_size * numbers_per_byte);
  else
  {
    // The default segment size is x^(1/4).
    // This is tiny, will fit into the CPU's cache.
    segment_size_ = x14_;
  }

  validate_segment_size();
}

bool LoadBalancerAC::get_work(int64_t& low, int64_t& high)
{
  LockGuard lockGuard(lock_);

  // Most special leaves are below y (~ x^(1/3) * log(x)).
  // We make sure this interval is evenly distributed
  // amongst all threads. Above y we slowly increase the
  // segment size but still ensure that it fits into the
  // CPU's cache.
  if (low_ > y_)
  {
    if (segment_size_ <= l2_cache_size * numbers_per_byte &&
        low_ + (l2_cache_size * numbers_per_byte * threads_) / 4 <= sqrtx_)
      segment_size_ = l2_cache_size * numbers_per_byte;
    else if (segment_size_ <= l1_cache_size * numbers_per_byte &&
              low_ + (l1_cache_size * numbers_per_byte * threads_) / 2 <= sqrtx_)
      segment_size_ = l1_cache_size * numbers_per_byte;
    else if (segment_size_ * 4 <= l1_cache_size * numbers_per_byte &&
              low_ + (segment_size_ * 4 * threads_) / 2 <= sqrtx_)
      segment_size_ *= 4;
  }

  validate_segment_size();
  low = low_;
  high = std::min(low + segment_size_, sqrtx_);
  low_ = high;

  return low < sqrtx_;
}

void LoadBalancerAC::validate_segment_size()
{
  segment_size_ = std::max(min_segment_size, segment_size_);

  if (segment_size_ % 240)
    segment_size_ += 240 - segment_size_ % 240;
}

} // namespace
