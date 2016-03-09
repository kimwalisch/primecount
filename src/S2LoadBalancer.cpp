///
/// @file  S2LoadBalancer.cpp
/// @brief The S2LoadBalancer evenly distributes the work load between
///        the threads in the computation of the special leaves.
///
/// Simply parallelizing the computation of the special leaves in the
/// Lagarias-Miller-Odlyzko and Deleglise-Rivat algorithms by
/// subdividing the sieve interval by the number of threads into
/// equally sized subintervals does not scale because the distribution
/// of the special leaves is highly skewed and most special leaves are
/// in the first few segments whereas later on there are very few
/// special leaves.
///
/// Based on the above observations it is clear that we need some kind
/// of load balancing in order to scale our parallel algorithm for
/// computing special leaves. Below are the ideas I used to develop a
/// load balancing algorithm that achieves a high load balance by
/// dynamically increasing or decreasing the interval size based on
/// the relative standard deviation of the thread run-times.
///
/// 1) Start with a tiny segment size of x^(1/3) / (log x * log log x)
///    and one segment per thread. Our algorithm uses equally sized
///    intervals, for each thread the interval_size is
///    segment_size * segments_per_thread and the threads process
///    adjacent intervals i.e.
///    [base + interval_size * thread_id, base + interval_size * (thread_id + 1)].
///
/// 2) If the relative standard deviation of the thread run-times is
///    large then we know the special leaves are distributed unevenly,
///    else if the relative standard deviation is low the special
///    leaves are more evenly distributed.
///
/// 3) If the special leaves are distributed unevenly then we can
///    increase the load balance by decreasing the interval_size.
///    Contrary if the special leaves are more evenly distributed
///    we can increase the interval_size in order to improve the
///    algorithm's efficiency.
///
/// 4) We can't use a static threshold for as to when the relative
///    standard deviation is low or large as this threshold varies for
///    different PC architectures. So instead we compare the current
///    relative standard deviation to the previous one in order to
///    decide whether to increase or decrease the interval_size.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
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

S2LoadBalancer::S2LoadBalancer(maxint_t x,
                               int64_t y,
                               int64_t z,
                               int64_t threads)
{
  double rsd = 40;
  init(x, y, z, threads, rsd);
}

S2LoadBalancer::S2LoadBalancer(maxint_t x,
                               int64_t y,
                               int64_t z,
                               int64_t threads,
                               double rsd)
{
  init(x, y, z, threads, rsd);
}

void S2LoadBalancer::init(maxint_t x,
                          int64_t y,
                          int64_t z,
                          int64_t threads,
                          double rsd)
{
  x_ = (double) x;
  y_ = (double) y;
  z_ = (double) z;
  rsd_ = rsd;
  avg_seconds_ = 0;
  count_ = 0;
  sqrtz_ = isqrt(z);

  // determined by benchmarking
  double log_threads = max(1.0, log((double) threads));
  decrease_dividend_ = max(0.5, log_threads / 3);
  min_seconds_ = 0.02 * log_threads;
  update_min_size(log(x_) * log(log(x_)));

  double alpha = get_alpha(x, y);
  smallest_hard_leaf_ = (int64_t) (x / (y * sqrt(alpha) * iroot<6>(x)));
}

double S2LoadBalancer::get_rsd() const
{
  return rsd_;
}

int64_t S2LoadBalancer::get_min_segment_size() const
{
  return min_size_;
}

bool S2LoadBalancer::decrease_size(double seconds,
                                   double decrease) const
{
  return seconds > min_seconds_ &&
         rsd_ > decrease;
}

bool S2LoadBalancer::increase_size(double seconds,
                                   double decrease) const
{
  return seconds < avg_seconds_ &&
        !decrease_size(seconds, decrease);
}

/// Used to decide whether to use a smaller or larger
/// segment_size and/or segments_per_thread.
///
double S2LoadBalancer::get_decrease_threshold(double seconds) const
{
  double log_seconds = max(min_seconds_, log(seconds));
  double dont_decrease = min(decrease_dividend_ / (seconds * log_seconds), rsd_);
  return rsd_ + dont_decrease;
}

void S2LoadBalancer::update_avg_seconds(double seconds)
{
  seconds = max(seconds, min_seconds_);
  double dividend = avg_seconds_ * count_ + seconds;
  avg_seconds_ = dividend / ++count_;
}

void S2LoadBalancer::update_min_size(double divisor)
{
  int64_t size = (int64_t) (sqrt(z_) / max(1.0, divisor));
  int64_t min_size = 1 << 9;
  min_size_ = max(size, min_size);
  min_size_ = next_power_of_2(min_size_);
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
  double decrease_threshold = get_decrease_threshold(seconds);
  rsd_ = max(0.1, relative_standard_deviation(timings));

  // 1 segment per thread
  if (*segment_size < sqrtz_)
  {
    if (increase_size(seconds, decrease_threshold))
      *segment_size <<= 1;
    else if (decrease_size(seconds, decrease_threshold))
      if (*segment_size > min_size_)
        *segment_size >>= 1;
  }
  // many segments per thread
  else if (low > smallest_hard_leaf_)
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

  int64_t high = low + *segment_size * *segments_per_thread * threads;

  // near smallest_hard_leaf_ the hard special leaves
  // are distributed unevenly so use min_size_
  if (low <= smallest_hard_leaf_ && 
      high > smallest_hard_leaf_)
  {
    *segment_size = min_size_;
  }

  high = low + *segment_size * *segments_per_thread * threads; 

  // slightly increase min_size_
  if (high >= smallest_hard_leaf_)
  {
    update_min_size(log(y_));
    *segment_size = max(min_size_, *segment_size);
  }
}

} // namespace
