///
/// @file  balance_S2_load.cpp
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

#include <primecount-internal.hpp>
#include <aligned_vector.hpp>
#include <pmath.hpp>

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
  double sum_mean_squared = 0;

  if (average == 0)
    return 0;

  for (size_t i = 0; i < n; i++)
  {
    double mean = timings[i] - average;
    sum_mean_squared += mean * mean;
  }

  double divisor = max(1.0, (double) n - 1);
  double standard_deviation = sqrt(sum_mean_squared / divisor);
  double rsd = 100 * standard_deviation / average;

  return rsd;
}

/// Used to decide whether to use a smaller or larger
/// segment_size and/or segments_per_thread.
///
double compute_decrease_threshold(double rsd, double seconds, double threads)
{
  double dividend = max(0.7, log(threads) / 3.0);
  double quotient = max(1.0, dividend / (seconds * log(seconds)));
  double dont_decrease = min(quotient, dividend * 10.0);
  return rsd + dont_decrease;
}

/// @param rsd  Relative standard deviation
bool decrease_size(double rsd,
                   double decrease_threshold,
                   double seconds)
{
  return seconds > 0.01 && rsd > decrease_threshold;
}

bool increase_size(double rsd,
                   double decrease_threshold,
                   double seconds,
                   double max_seconds)
{
  return seconds < max_seconds && (seconds < 0.01 || rsd < decrease_threshold);
}

bool adjust_segments(double segments,
                     double seconds,
                     double max_seconds,
                     int64_t old_segments)
{
  return ((segments < old_segments) && (seconds > 0.01)) ||
         ((segments > old_segments) && (seconds < max_seconds));
}

} // namespace

namespace primecount {

/// Balance the load in the computation of the special leaves
/// by dynamically adjusting the segment_size and segments_per_thread.
/// @param old_rsd  Previous relative standard deviation.
/// @param timings  Timings of the threads.
///
void balance_S2_load(double x,
                     double threads,
                     double* old_rsd,
                     aligned_vector<double>& timings,
                     int64_t* segment_size,
                     int64_t* segments_per_thread,
                     int64_t min_segment_size,
                     int64_t max_segment_size)
{
  double seconds = get_average(timings);
  double rsd = max(0.1, relative_standard_deviation(timings));
  double max_seconds = in_between(5, log10(x) * max(1.0, log10(threads)), 60);
  double decrease_threshold = compute_decrease_threshold(*old_rsd, seconds, threads);

  if (*segment_size < max_segment_size)
  {
    if (increase_size(rsd, decrease_threshold, seconds, max_seconds))
      *segment_size <<= 1;
    else if (decrease_size(rsd, decrease_threshold, seconds))
      if (*segment_size > min_segment_size)
        *segment_size >>= 1;
  }
  else
  {
    double factor = decrease_threshold / rsd;
    factor = in_between(0.25, factor, 4.0);
    double segments = *segments_per_thread * factor;
    segments = max(1.0, segments);

    if (adjust_segments(segments, seconds, max_seconds, *segments_per_thread))
      *segments_per_thread = (int) segments;
  }

  *old_rsd = rsd;
}

} //namespace
