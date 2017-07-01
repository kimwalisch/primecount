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
#include <Sieve.hpp>
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
                           maxint_t s2_approx) :
  low_(1),
  max_low_(1),
  z_(z),
  segments_(1),
  s2_total_(0),
  s2_approx_(s2_approx),
  time_(get_wtime()),
  status_(x)
{
  init_size();
  maxint_t x16 = iroot<6>(x);
  double alpha = get_alpha(x, y);
  smallest_hard_leaf_ = (int64_t) (x / (y * sqrt(alpha) * x16));
}

void LoadBalancer::init_size()
{
  // start with a tiny segment_size as most
  // special leaves are in the first few segments
  // and we need to ensure that all threads are
  // assigned an equal amount of work
  int64_t sqrtz = isqrt(z_);
  int64_t log = ilog(sqrtz);
  log = max(log, 1);
  segment_size_ = sqrtz / log;

  int64_t min_size = 1 << 9;
  segment_size_ = max(segment_size_, min_size);
  segment_size_ = Sieve::get_segment_size(segment_size_);

  // try to use a segment size that fits exactly
  // into the CPUs L1 data cache
  int64_t l1_dcache_size = 1 << 15;
  max_size_ = l1_dcache_size * 30;
  max_size_ = max(max_size_, sqrtz);
  max_size_ = Sieve::get_segment_size(max_size_);
}

maxint_t LoadBalancer::get_result() const
{
  return s2_total_;
}

bool LoadBalancer::get_work(int64_t* low,
                            int64_t* segments,
                            int64_t* segment_size,
                            maxint_t s2,
                            Runtime& runtime)
{
  #pragma omp critical (get_work)
  {
    s2_total_ += s2;

    update(low, segments, runtime);

    if (low_ == 1)
      segment_size_ -= 1;
    else
      segment_size_ = Sieve::get_segment_size(segment_size_);

    *low = low_;
    *segments = segments_;
    *segment_size = segment_size_;

    low_ += segments_ * segment_size_;

    if (is_print())
      status_.print(s2_total_, s2_approx_);
  }

  return *low <= z_;
}

void LoadBalancer::update(int64_t* low,
                          int64_t* segments,
                          Runtime& runtime)
{
  if (*low > max_low_)
  {
    max_low_ = *low;
    segments_ = *segments;

    if (segment_size_ < max_size_)
      segment_size_ = min(segment_size_ * 2, max_size_);
    else
    {
      double next = get_next(runtime);
      next = in_between(0.25, next, 2.0);
      next *= segments_;
      next = max(1.0, next);
      segments_ = (int64_t) next;
    }
  }

  auto high = low_ + segments_ * segment_size_;

  // Most hard special leaves are located just past
  // smallest_hard_leaf_. In order to prevent assigning
  // the bulk of work to a single thread we reduce
  // the number of segments to a minimum.
  //
  if (smallest_hard_leaf_ >= low_ &&
      smallest_hard_leaf_ <= high)
  {
    segments_ = 1;
  }
}

double LoadBalancer::get_next(Runtime& runtime) const
{
  double min_secs = runtime.init * 10;
  double run_secs = runtime.secs;

  min_secs = max(min_secs, 0.01);
  run_secs = max(run_secs, min_secs / 10);

  double rem = remaining_secs();
  double threshold = rem / 4;
  threshold = max(threshold, min_secs);

  return threshold / run_secs;
}

/// Remaining seconds till finished
double LoadBalancer::remaining_secs() const
{
  double percent = status_.getPercent(low_, z_, s2_total_, s2_approx_);
  percent = in_between(20, percent, 100);

  double total_secs = get_wtime() - time_;
  double secs = total_secs * (100 / percent) - total_secs;
  return secs;
}

} // namespace
