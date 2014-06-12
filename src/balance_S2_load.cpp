///
/// @file  balance_S2_load.cpp
/// @brief Dynamically increase or decrease the segment_size or the
///        segments_per_thread in order to improve the load balancing
///        in the computation of the special leaves. The algorithm
///        calculates the relative standard deviation of the timings
///        of the individual threads and then decides whether to
///        assign more or less work to the threads.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <balance_S2_load.hpp>
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
    // error correction
    mean = min(mean, average * 2);
    sum_mean_squared += mean * mean;
  }

  double standard_deviation = sqrt(sum_mean_squared / n);
  double rel_standard_deviation = 100 * standard_deviation / average;

  return rel_standard_deviation;
}

} // namespace

namespace primecount {

void balance_S2_load(int64_t* segment_size,
                     int64_t* segments_per_thread,
                     int64_t min_segment_size,
                     int64_t max_segment_size,
                     aligned_vector<double>& timings)
{
  double seconds = get_average(timings);
  double rsd = relative_standard_deviation(timings);
  double increase_threshold = 6;

  if (seconds < 0.1)
    increase_threshold *= 2;
  else if (seconds > 2)
    increase_threshold /= 2;

  rsd = in_between(increase_threshold / 7, rsd, increase_threshold * 7);

  if (*segment_size < max_segment_size)
  {
    // balance segment_size
    if ((seconds < 0.1 || rsd < increase_threshold) &&
         seconds < 10)
      *segment_size <<= 1;
  }
  else
  {
    // balance segments_per_thread
    double segments = max(1.0, *segments_per_thread * increase_threshold / rsd);

    if ((segments < *segments_per_thread && seconds > 0.01) ||
        (segments > *segments_per_thread && seconds < 10))
    {
      *segments_per_thread = (int) segments;
    }
  }
}

} //namespace
