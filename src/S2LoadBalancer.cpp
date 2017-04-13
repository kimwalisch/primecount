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
///    and one segment per thread. Each thread sieves a distance of
///    segment_size * segments_per_thread. Once all threads have
///    processed their intervals we calculate a new interval size
///    based on the ideas below.
///
/// 2) If the relative standard deviation of the thread run-times is
///    large then we know the special leaves are distributed unevenly,
///    else if the relative standard deviation is low the special
///    leaves are more evenly distributed.
///
/// 3) If the special leaves are distributed unevenly then we can
///    increase the load balance by decreasing the interval size.
///    Contrary if the special leaves are more evenly distributed
///    we can increase the interval size in order to improve the
///    algorithm's efficiency.
///
/// 4) We can't use a static threshold for as to when the relative
///    standard deviation is low or large as this threshold varies for
///    different PC architectures. So instead we compare the current
///    relative standard deviation to the previous one in order to
///    decide whether to increase or decrease the interval size.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <S2LoadBalancer.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <algorithm>
#include <cmath>

using namespace std;
using namespace primecount;

namespace {

double get_avg(thread_timings_t& timings)
{
  size_t n = timings.size();
  double sum = 0;

  for (size_t i = 0; i < n; i++)
    sum += timings[i];

  return sum / n;
}

/// Relative standard deviation
double rel_std_dev(thread_timings_t& timings)
{
  size_t n = timings.size();
  double avg = get_avg(timings);
  double sum = 0;

  if (avg == 0)
    return 0;

  for (size_t i = 0; i < n; i++)
  {
    double mean = timings[i] - avg;
    sum += mean * mean;
  }

  double std_dev = sqrt(sum / max(1.0, n - 1.0));
  double rsd = 100 * std_dev / avg;

  return rsd;
}

} // namespace

namespace primecount {

S2LoadBalancer::S2LoadBalancer(maxint_t x,
                               int64_t y,
                               int64_t z,
                               int64_t threads) :
  x_((double) x),
  y_((double) y),
  z_((double) z),
  rsd_(40),
  count_(0),
  total_seconds_(0),
  sqrtz_(isqrt(z))
{
  init(x, y, threads);
}

S2LoadBalancer::S2LoadBalancer(maxint_t x,
                               int64_t y,
                               int64_t z,
                               int64_t threads,
                               double rsd) :
  x_((double) x),
  y_((double) y),
  z_((double) z),
  rsd_(rsd),
  count_(0),
  total_seconds_(0),
  sqrtz_(isqrt(z))
{
  init(x, y, threads);
}

void S2LoadBalancer::init(maxint_t x,
                          int64_t y,
                          int64_t threads)
{
  // determined by benchmarking
  double log_threads = max(1.0, log((double) threads));
  decrease_dividend_ = max(0.5, log_threads / 3);

  min_seconds_ = 0.01 * log_threads;
  double divisor = log(log(x_)) * log(x_);
  update_min_size(divisor);

  double alpha = get_alpha(x, y);
  smallest_hard_leaf_ = (int64_t) (x / (y * sqrt(alpha) * iroot<6>(x)));
}

double S2LoadBalancer::get_rsd() const
{
  return rsd_;
}

double S2LoadBalancer::get_avg_seconds() const
{
  return total_seconds_ / count_;
}

int64_t S2LoadBalancer::get_min_segment_size() const
{
  return min_size_;
}

/// Increase if relative std dev < pivot.
/// Decrease if relative std dev > pivot.
///
double S2LoadBalancer::get_pivot(double seconds) const
{
  double log_seconds = log(seconds);
  log_seconds = max(min_seconds_, log_seconds);
  double dont_decrease = decrease_dividend_ / (seconds * log_seconds);
  dont_decrease = min(dont_decrease, rsd_);

  return rsd_ + dont_decrease;
}

bool S2LoadBalancer::is_increase(double seconds,
                                 double pivot) const
{
  return (seconds < min_seconds_ ||
          seconds < get_avg_seconds()) &&
         !is_decrease(seconds, pivot);
}

bool S2LoadBalancer::is_decrease(double seconds,
                                 double pivot) const
{
  return seconds > min_seconds_ &&
         rsd_ > pivot;
}

void S2LoadBalancer::update_min_size(double divisor)
{
  int64_t min_size = 1 << 9;
  int64_t size = (int64_t) (sqrtz_ / max(1.0, divisor));
  min_size_ = max(size, min_size);
  min_size_ = next_power_of_2(min_size_);
}

void S2LoadBalancer::update(int64_t* segment_size,
                            int64_t* segments_per_thread,
                            int64_t low,
                            int64_t threads,
                            thread_timings_t& timings)
{
  count_++;
  double seconds = get_avg(timings);
  total_seconds_ += seconds;
  double pivot = get_pivot(seconds);
  rsd_ = max(0.1, rel_std_dev(timings));

  // 1 segment per thread
  if (*segment_size < sqrtz_)
  {
    if (is_increase(seconds, pivot))
      *segment_size <<= 1;
    else if (is_decrease(seconds, pivot))
      if (*segment_size > min_size_)
        *segment_size >>= 1;
  }
  // many segments per thread
  else if (low > smallest_hard_leaf_)
    update(segments_per_thread, seconds, pivot);

  int64_t thread_distance = *segment_size * *segments_per_thread;
  int64_t high = low + thread_distance * threads;

  // near smallest_hard_leaf_ the hard special leaves
  // are distributed unevenly so use min_size_
  if (low <= smallest_hard_leaf_ && 
      high > smallest_hard_leaf_)
  {
    *segment_size = min_size_;
    thread_distance = *segment_size * *segments_per_thread;
    high = low + thread_distance * threads;
  }

  // slightly increase min_size_
  if (high >= smallest_hard_leaf_)
  {
    update_min_size(log(y_));
    *segment_size = max(min_size_, *segment_size);
  }
}

/// Increase the segments_per_thread if the relative standard
/// deviation of the thread run times is small, or
/// decrease the segments_per_thread if the relative standard
/// deviation of the thread run times is large.
///
void S2LoadBalancer::update(int64_t* segments_per_thread,
                            double seconds,
                            double pivot)
{
  if (is_increase(seconds, pivot) ||
      is_decrease(seconds, pivot))
  {
    if (seconds < min_seconds_)
      *segments_per_thread *= 2;
    else
    {
      double factor = pivot / rsd_;
      factor = in_between(0.5, factor, 2);
      double segments = *segments_per_thread * factor;
      segments = max(1.0, segments);
      *segments_per_thread = (int64_t) segments;
    }
  }
}

} // namespace
