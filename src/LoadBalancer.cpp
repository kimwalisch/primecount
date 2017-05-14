///
/// @file  LoadBalancer.cpp
/// @brief The LoadBalancer assigns work to the individual
///        threads in the computation of the special leaves
///        in the Lagarias-Miller-Odlyzko and
///        Deleglise-Rivat prime counting algorithms.
///
/// Simply parallelizing the computation of the special leaves in the
/// Lagarias-Miller-Odlyzko and Deleglise-Rivat algorithms by
/// subdividing the sieve interval by the number of threads into
/// equally sized subintervals does not scale because the distribution
/// of the special leaves is highly skewed and most special leaves are
/// in the first few segments whereas later on there are very few
/// special leaves.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <LoadBalancer.hpp>
#include <primecount-internal.hpp>
#include <S2Status.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <min.hpp>

#include <stdint.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

namespace primecount {

LoadBalancer::LoadBalancer(maxint_t x, int64_t y, maxint_t s2_approx) :
  status_(x),
  low_(1),
  limit_(x / y + 1),
  segments_(1),
  s2_approx_(s2_approx),
  S2_total_(0),
  time_(get_wtime()),
  finished_(false)
{
  init_size(x, y);
}

void LoadBalancer::init_size(maxint_t x, int64_t y)
{
  double n = (double) x;
  double divisor = log(log(n)) * log(n);

  int64_t min_size = 1 << 9;
  int64_t z = (int64_t) (x / y);
  int64_t sqrtz = isqrt(z);

  // Start with a tiny segment size
  segment_size_ = (int64_t) (sqrtz / max(1.0, divisor));
  segment_size_ = max(min_size, segment_size_);
  segment_size_ = next_power_of_2(segment_size_);

  max_size_ = max(min_size, sqrtz);
  max_size_ = next_power_of_2(max_size_);
}

bool LoadBalancer::is_increase(double percent, Runtime& runtime)
{
  double min_seconds = max(0.01, runtime.init * 10);

  if (runtime.seconds < min_seconds)
    return true;

  // avoid division by 0
  percent = in_between(1, percent, 99.9);

  // remaining time till finished
  double total_time = get_wtime() - time_;
  double remaining_time = total_time * (100 / percent) - total_time;
  double threshold = remaining_time / 4;
  threshold = max(min_seconds, threshold);

  return runtime.seconds < threshold;
}

void LoadBalancer::get_work(int64_t* low,
                            int64_t* segments,
                            int64_t* segment_size,
                            maxint_t S2, 
                            Runtime& runtime)
{
  #pragma omp critical (S2_schedule)
  {
    *low = low_;
    *segments = segments_;
    *segment_size = segment_size_;

    S2_total_ += S2;
    low_ += segments_ * segment_size_;

    if (*low >= limit_)
      finished_ = true;
    else
    {
      double percent = status_.skewed_percent(S2_total_, s2_approx_);

      if (!is_increase(percent, runtime))
        segments_ -= segments_ / 4;
      else
      {
        if (segment_size_ < max_size_)
          segment_size_ *= 2;
        else
          segments_ += max(segments_ / 3, 1);
      }
    }
  }

  if (is_print())
    status_.print(S2_total_, s2_approx_);
}

bool LoadBalancer::finished()
{
  return finished_;
}

maxint_t LoadBalancer::get_result()
{
  return S2_total_;
}

} // namespace
