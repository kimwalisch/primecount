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
#include <cmath>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

namespace primecount {

LoadBalancer::LoadBalancer(maxint_t x,
                          int64_t y,
                          int64_t z,
                          double alpha,
                          maxint_t s2_approx) :
  low_(1),
  max_low_(1),
  z_(z),
  segments_(1),
  s2_approx_(s2_approx),
  S2_total_(0),
  time_(get_wtime()),
  status_(x)
{
  init_size();
  maxint_t x16 = iroot<6>(x);
  smallest_hard_leaf_ = (int64_t) (x / (y * sqrt(alpha) * x16));
}

void LoadBalancer::init_size()
{
  int64_t min_size = 1 << 9;
  int64_t sqrtz = isqrt(z_);
  segment_size_ = max(min_size, sqrtz);
  segment_size_ = next_power_of_2(segment_size_);
}

maxint_t LoadBalancer::get_result() const
{
  return S2_total_;
}

bool LoadBalancer::get_work(int64_t* low,
                            int64_t* segments,
                            int64_t* segment_size,
                            maxint_t S2,
                            Runtime& runtime)
{
  #pragma omp critical (S2_schedule)
  {
    update(low, segments, segment_size, runtime);

    *low = low_;
    *segments = segments_;
    *segment_size = segment_size_;

    S2_total_ += S2;
    low_ += segments_ * segment_size_;

    if (is_print())
      status_.print(max_low_, z_, S2_total_, s2_approx_);
  }

  return *low <= z_;
}

void LoadBalancer::update(int64_t* low,
                          int64_t* segments,
                          int64_t* segment_size,
                          Runtime& runtime)
{
  if (*low > max_low_)
  {
    max_low_ = *low;
    segments_ = *segments;
    segment_size_ = *segment_size;

    if (is_increase(runtime))
      segments_ *= 2;
    else
      segments_ = ceil_div(segments_, 2);
  }

  // Most hard special leaves are located just past
  // smallest_hard_leaf_. In order to prevent assigning
  // the bulk of work to a single thread we reduce
  // the number of segments to a minimum.
  //
  if (smallest_hard_leaf_ >= low_ &&
      smallest_hard_leaf_ <= low_ + segments_ * segment_size_)
  {
    segments_ = 1;
  }
}

bool LoadBalancer::is_increase(Runtime& runtime) const
{
  double min_secs = max(0.01, runtime.init * 10);

  if (runtime.secs < min_secs)
    return true;

  double secs = remaining_secs();
  double threshold = secs / 4;
  threshold = max(min_secs, threshold);

  return runtime.secs < threshold;
}

/// Remaining seconds till finished
double LoadBalancer::remaining_secs() const
{
  double percent = status_.getPercent(low_, z_, S2_total_, s2_approx_);
  percent = in_between(10, percent, 99.9);

  double total_secs = get_wtime() - time_;
  double secs = total_secs * (100 / percent) - total_secs;
  return secs;
}

} // namespace
