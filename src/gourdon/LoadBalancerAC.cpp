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
#include <imath.hpp>

#include <stdint.h>
#include <cmath>

using namespace std;

namespace primecount {

LoadBalancerAC::LoadBalancerAC(int64_t sqrtx,
                               int64_t x13,
                               int64_t y,
                               int threads) :
  low_(0),
  sqrtx_(sqrtx),
  x13_(x13),
  x14_(isqrt(sqrtx)),
  x29_(std::pow(x13, 2.0 / 3)),
  y_(y),
  threads_(threads)
{ }

bool LoadBalancerAC::get_work(int64_t& low, int64_t& high)
{
  LockGuard lockGuard(lock_);

  // CPU cache sizes per core
  int64_t l1_cache_size = 32 << 10;
  int64_t l2_cache_size = 256 << 10;

  // numbers_per_byte = 240 / sizeof(SegmentedPiTable::pi_t)
  int64_t numbers_per_byte = 15;
  int64_t segment_size;

  if (threads_ == 1)
    segment_size = max(x14_, l2_cache_size * numbers_per_byte);
  else
  {
    // The default segment size is x^(1/4).
    // This is tiny, will fit into the CPU's cache.
    segment_size = x14_;

    // Most special leaves are below x^(1/3).
    // We make sure this interval is evenly distributed
    // amongst all threads. Above x^(1/3) we slowly
    // increase the segment size but still ensure that
    // it fits into the CPU's cache.
    if (low_ <= x13_ && x14_ * threads_ > x13_)
      segment_size = x29_;
    else if (low_ <= x13_ * 4)
      segment_size = x14_;
    else if (low_ <= y_)
      segment_size = x14_ * 2;
    else if (segment_size / numbers_per_byte <= l2_cache_size &&
              low_ + (l2_cache_size * numbers_per_byte * threads_) / 4 <= sqrtx_)
      segment_size = l2_cache_size * numbers_per_byte;
    else if (segment_size / numbers_per_byte <= l1_cache_size &&
              low_ + (l1_cache_size * numbers_per_byte * threads_) / 2 <= sqrtx_)
      segment_size = l1_cache_size * numbers_per_byte;
  }

  // Minimum segment size = 1 KiB
  int64_t min_segment_size = (1 << 10) * numbers_per_byte;
  segment_size = max(min_segment_size, segment_size);

  if (segment_size % 240)
    segment_size += 240 - segment_size % 240;

  low = low_;
  high = std::min(low + segment_size, sqrtx_);
  low_ += segment_size;

  return low < sqrtx_;
}

} // namespace
