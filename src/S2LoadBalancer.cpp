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
#include <ptypes.hpp>

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

S2LoadBalancer::S2LoadBalancer(maxint_t x, int64_t z) :
  x_((double) x),
  rsd_(40)
{
  double divisor = log(x_) * log(log(x_));
  int64_t size = (int64_t) (isqrt(z) / max(1.0, divisor));
  size = max((int64_t) (1 << 9), size);

  min_size_ = next_power_of_2(size);
  max_size_ = next_power_of_2(isqrt(z));
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
  return seconds > 0.01 && rsd_ > decrease;
}

bool S2LoadBalancer::increase_size(double seconds, double max_seconds, double decrease) const
{
  return seconds < max_seconds &&
        !decrease_size(seconds, decrease);
}

/// Synchronize threads after at most max_seconds
double S2LoadBalancer::get_max_seconds(int64_t threads) const
{
  double log_threads = log10((double) threads);
  return max(2.0, log10(x_) * max(1.0, log_threads));
}

/// Used to decide whether to use a smaller or larger
/// segment_size and/or segments_per_thread.
///
double S2LoadBalancer::get_decrease_threshold(double seconds, int64_t threads) const
{
  double t = (double) threads;
  double dividend = max(0.5, log(t) / 4);
  double quotient = max(dividend, dividend / (seconds * log(seconds)));
  double dont_decrease = min(quotient, dividend * 10);
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
  double max_seconds = get_max_seconds(threads);
  double decrease_threshold = get_decrease_threshold(seconds, threads);
  rsd_ = max(0.1, relative_standard_deviation(timings));

  // 1 segment per thread
  if (*segment_size < max_size_)
  {
    if (increase_size(seconds, max_seconds, decrease_threshold))
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
    factor = in_between(0.25, factor, 4.0);
    double n = *segments_per_thread * factor;
    n = max(1.0, n);

    if ((n < *segments_per_thread && seconds > 0.01) ||
        (n > *segments_per_thread && seconds < max_seconds))
      *segments_per_thread = (int) n;
  }
}

} //namespace
