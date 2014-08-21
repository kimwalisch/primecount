///
/// @file  S2LoadBalancer.cpp
/// @brief Dynamically increase or decrease the segment_size or the
///        segments_per_thread in order to improve the load balancing
///        in the computation of the special leaves. The algorithm
///        calculates the relative standard deviation of the timings
///        of the individual threads and then decides whether to
///        assign more work (low relative standard deviation) or less
///        work (large relative standard deviation) to the threads.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <S2LoadBalancer.hpp>
#include <primecount-internal.hpp>
#include <aligned_vector.hpp>
#include <pmath.hpp>
#include <int128.hpp>

#include <stdint.h>
#include <algorithm>
#include <cmath>

using namespace std;
using namespace primecount;

namespace {

double get_average(aligned_vector<double>& timings)
{
  size_t n = timings.size();
  double sum = 0;

  for (size_t i = 0; i < n; i++)
    sum += timings[i];

  return sum / n;
}

double relative_standard_deviation(aligned_vector<double>& timings)
{
  size_t n = timings.size();
  double average = get_average(timings);
  double sum = 0;

  if (average == 0)
    return 0;

  for (size_t i = 0; i < n; i++)
  {
    double mean = timings[i] - average;
    sum += mean * mean;
  }

  double std_dev = sqrt(sum / max(1.0, n - 1.0));
  double rsd = 100 * std_dev / average;

  return rsd;
}

} // namespace

namespace primecount {

S2LoadBalancer::S2LoadBalancer(maxint_t x, int64_t z, int64_t threads) :
  x_((double) x),
  z_((double) z),
  rsd_(40),
  avg_seconds_(0),
  count_(0)
{
  double log_threads = max(1.0, log((double) threads));
  min_seconds_ = 0.02 * log_threads;
  max_size_ = next_power_of_2(isqrt(z));
  update_min_size(log(x_) * log(log(x_)));
}

double S2LoadBalancer::get_rsd() const
{
  return rsd_;
}

int64_t S2LoadBalancer::get_min_segment_size() const
{
  return min_size_;
}

bool S2LoadBalancer::decrease_size(double seconds, double decrease) const
{
  return seconds > min_seconds_ &&
         rsd_ > decrease;
}

bool S2LoadBalancer::increase_size(double seconds, double decrease) const
{
  return seconds < avg_seconds_ &&
        !decrease_size(seconds, decrease);
}

void S2LoadBalancer::update_avg_seconds(double seconds)
{
  seconds = max(seconds, min_seconds_);
  double dividend = avg_seconds_ * count_ + seconds;
  avg_seconds_ = dividend / (count_ + 1);
  count_++;
}

void S2LoadBalancer::update_min_size(double divisor)
{
  double size = sqrt(z_) / max(1.0, divisor);
  min_size_ = max((int64_t) (1 << 9), (int64_t) size);
  min_size_ = next_power_of_2(min_size_);
}

/// Used to decide whether to use a smaller or larger
/// segment_size and/or segments_per_thread.
///
double S2LoadBalancer::get_decrease_threshold(double seconds, int64_t threads) const
{
  double log_threads = log((double) threads);
  double log_seconds = max(min_seconds_, log(seconds));
  double dividend = max(0.1, log_threads / 4.0);
  double dont_decrease = min(dividend / (seconds * log_seconds), rsd_);
  return rsd_ + dont_decrease;
}

/// Balance the load in the computation of the special leaves
/// by dynamically adjusting the segment_size and segments_per_thread.
/// @param timings  Timings of the threads.
///
void S2LoadBalancer::update(int64_t low,
                            int64_t threads,
                            int64_t* segment_size,
                            int64_t* segments_per_thread,
                            aligned_vector<double>& timings)
{
  double seconds = get_average(timings);
  update_avg_seconds(seconds);
  double decrease_threshold = get_decrease_threshold(seconds, threads);
  rsd_ = max(0.1, relative_standard_deviation(timings));

  // if low > sqrt(z) we use a larger min_size_ as the
  // special leaves are distributed more evenly
  if (low > max_size_)
    update_min_size(log(x_));

  // 1 segment per thread
  if (*segment_size < max_size_)
  {
    if (increase_size(seconds, decrease_threshold))
      *segment_size <<= 1;
    else if (decrease_size(seconds, decrease_threshold))
      if (*segment_size > min_size_)
        *segment_size >>= 1;

    // near sqrt(z) there is a short peak of special
    // leaves so we use the minium segment size
    int64_t high = low + *segment_size * *segments_per_thread * threads; 
    if (low <= max_size_ && high > max_size_)
      *segment_size = min_size_;
  }
  else // many segments per thread
  {
    double factor = decrease_threshold / rsd_;
    factor = in_between(0.5, factor, 2);
    double n = *segments_per_thread * factor;
    n = max(1.0, n);

    if ((n < *segments_per_thread && seconds > min_seconds_) ||
        (n > *segments_per_thread && seconds < avg_seconds_))
    {
      *segments_per_thread = (int) n;
    }
  }
}

} //namespace
