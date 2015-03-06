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
/// Copyright (C) 2015 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <S2LoadBalancer.hpp>
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
  decrease_dividend_ = max(0.5, log_threads / 3);
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
  avg_seconds_ = dividend / (count_ + 1);
  count_++;
}

void S2LoadBalancer::update_min_size(double divisor)
{
  double size = sqrt(z_) / max(1.0, divisor);
  min_size_ = max((int64_t) (1 << 9), (int64_t) size);
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

  // if low > sqrt(z) we use a larger min_size_ as the
  // special leaves are distributed more evenly
  if (low > max_size_)
  {
    update_min_size(log(x_));
    *segment_size = max(*segment_size, min_size_);
  }

  // 1 segment per thread
  if (*segment_size < max_size_)
  {
    if (increase_size(seconds, decrease_threshold))
      *segment_size <<= 1;
    else if (decrease_size(seconds, decrease_threshold))
      if (*segment_size > min_size_)
        *segment_size >>= 1;

    // near sqrt(z) there is a short peak of special
    // leaves so we use the minimum segment size
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

} // namespace
