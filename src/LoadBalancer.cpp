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
#include <calculator.hpp>
#include <json.hpp>

#include <stdint.h>
#include <cmath>
#include <iostream>
#include <iomanip>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace nlohmann;

namespace primecount {

int LoadBalancer::resume_threads() const
{
  int threads = 0;
  auto j = load_backup();

  if (is_resume(j, "S2_hard", x_, y_, z_))
    while (j["S2_hard"]["thread" + to_string(threads)].is_object())
      threads++;

  return threads;
}

void LoadBalancer::backup(int thread_id,
                          int64_t low,
                          int64_t segments,
                          int64_t segment_size) const
{
  auto j = load_backup();

  double percent = status_.getPercent(low_, z_, s2_total_, s2_approx_);

  j["S2_hard"]["x"] = to_string(x_);
  j["S2_hard"]["y"] = y_;
  j["S2_hard"]["z"] = z_;
  j["S2_hard"]["low"] = low_;
  j["S2_hard"]["segments"] = segments_;
  j["S2_hard"]["segment_size"] = segment_size_;
  j["S2_hard"]["s2_hard"] = to_string(s2_total_);
  j["S2_hard"]["percent"] = percent;
  j["S2_hard"]["seconds"] = get_wtime() - time_;

  j["S2_hard"]["thread" + to_string(thread_id)]["low"] = low;
  j["S2_hard"]["thread" + to_string(thread_id)]["segments"] = segments;
  j["S2_hard"]["thread" + to_string(thread_id)]["segment_size"] = segment_size;

  store_backup(j);
}

template <typename T>
void print_resume(int thread_id,
                  T low,
                  T segments,
                  T segment_size)
{
  if (is_print())
  {
    if (!print_variables())
      cout << endl;

    cout << "=== Resuming from primecount.backup ===" << endl;
    cout << "thread = " << thread_id << endl;
    cout << "low = " << low << endl;
    cout << "segments = " << segments << endl;
    cout << "segment_size = " << segment_size << endl;
    cout << endl;
  }
}

template <typename T>
void print_resume(T s2_hard,
                  double seconds)
{
  if (is_print())
  {
    if (!print_variables())
      cout << endl;

    cout << "=== Resuming from primecount.backup ===" << endl;
    cout << "S2_hard = " << s2_hard << endl;
    cout << "Seconds: " << fixed << setprecision(3) << seconds << endl;
    cout << endl;
  }
}

bool LoadBalancer::resume(int thread_id,
                          int64_t& low,
                          int64_t& segments,
                          int64_t& segment_size)
{
  bool resumed = false;

  #pragma omp critical (LoadBalancer)
  {
    auto j = load_backup();

    if (is_resume(j, "S2_hard", x_, y_, z_))
    {
      double seconds = j["S2_hard"]["seconds"];
      s2_total_ = calculator::eval<maxint_t>(j["S2_hard"]["s2_hard"]);
      low_ = j["S2_hard"]["low"];
      segments_ = j["S2_hard"]["segments"];
      segment_size_ = j["S2_hard"]["segment_size"];
      time_ = get_wtime() - seconds;
      low = j["S2_hard"]["thread" + to_string(thread_id)]["low"];
      segments = j["S2_hard"]["thread" + to_string(thread_id)]["segments"];
      segment_size = j["S2_hard"]["thread" + to_string(thread_id)]["segment_size"];

      print_resume(thread_id, low, segments, segment_size);

      resumed = true;
    }
  }

  return resumed;
}

bool LoadBalancer::resume(maxint_t& s2_hard, double& time) const
{
  auto j = load_backup();

  if (is_resume(j, "S2_hard", x_, y_, z_) &&
      j["S2_hard"]["low"] >= j["S2_hard"]["z"])
  {
    double seconds = j["S2_hard"]["seconds"];
    s2_hard = calculator::eval<maxint_t>(j["S2_hard"]["s2_hard"]);
    time = get_wtime() - seconds;
    print_resume(s2_hard, seconds);

    return true;
  }

  return false;
}

void LoadBalancer::finish_resume(int thread_id,
                                 int64_t low,
                                 int64_t segments,
                                 int64_t segment_size,
                                 maxint_t s2,
                                 Runtime& runtime)
{
  #pragma omp critical (LoadBalancer)
  {
    s2_total_ += s2;

    update(&low, &segments, runtime);

    double percent = status_.getPercent(low_, z_, s2_total_, s2_approx_);
    auto j = load_backup();

    j["S2_hard"]["x"] = to_string(x_);
    j["S2_hard"]["y"] = y_;
    j["S2_hard"]["z"] = z_;
    j["S2_hard"]["low"] = low_;
    j["S2_hard"]["segments"] = segments_;
    j["S2_hard"]["segment_size"] = segment_size_;
    j["S2_hard"]["s2_hard"] = to_string(s2_total_);
    j["S2_hard"]["percent"] = percent;
    j["S2_hard"]["seconds"] = get_wtime() - time_;
    j["S2_hard"].erase("thread" + to_string(thread_id));

    store_backup(j);

    if (is_print())
      status_.print(s2_total_, s2_approx_);
  }
}

void LoadBalancer::backup_result() const
{
  auto j = load_backup();

  if (is_resume(j, "S2_hard", x_, y_, z_))
  {
    j.erase("S2_hard");

    j["S2_hard"]["x"] = to_string(x_);
    j["S2_hard"]["y"] = y_;
    j["S2_hard"]["z"] = z_;
    j["S2_hard"]["low"] = low_;
    j["S2_hard"]["s2_hard"] = to_string(s2_total_);
    j["S2_hard"]["percent"] = 100;
    j["S2_hard"]["seconds"] = get_wtime() - time_;

    store_backup(j);
  }
}

LoadBalancer::LoadBalancer(maxint_t x,
                           int64_t y,
                           int64_t z,
                           double alpha,
                           maxint_t s2_approx) :
  low_(1),
  max_low_(1),
  x_(x),
  y_(y),
  z_(z),
  sqrtz_(isqrt(z)),
  segments_(1),
  s2_total_(0),
  s2_approx_(s2_approx),
  time_(get_wtime()),
  status_(x)
{
  init_size();
  maxint_t x16 = iroot<6>(x);
  smallest_hard_leaf_ = (int64_t) (x / (y * sqrt(alpha) * x16));
}

void LoadBalancer::init_size()
{
  int64_t log = ilog(sqrtz_);
  log = max(log, 1);
  segment_size_ = sqrtz_ / log;

  int64_t min = 1 << 9;
  segment_size_ = max(segment_size_, min);
  segment_size_ = next_power_of_2(segment_size_);
}

maxint_t LoadBalancer::get_result() const
{
  return s2_total_;
}

double LoadBalancer::get_time() const
{
  return time_;
}

bool LoadBalancer::get_work(int thread_id,
                            int64_t* low,
                            int64_t* segments,
                            int64_t* segment_size,
                            maxint_t s2,
                            Runtime& runtime)
{
  #pragma omp critical (LoadBalancer)
  {
    s2_total_ += s2;

    update(low, segments, runtime);

    *low = low_;
    *segments = segments_;
    *segment_size = segment_size_;

    low_ += segments_ * segment_size_;

    backup(thread_id, *low, *segments, *segment_size);

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

    if (segment_size_ < sqrtz_)
      segment_size_ *= 2;
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
