///
/// @file  LoadBalancer.cpp
/// @brief The LoadBalancer assigns work to the individual threads
///        in the computation of the special leaves in the
///        Lagarias-Miller-Odlyzko, Deleglise-Rivat and Gourdon
///        prime counting algorithms.
///
///        Simply parallelizing the computation of the special
///        leaves in the Lagarias-Miller-Odlyzko algorithm by
///        subdividing the sieve interval by the number of threads
///        into equally sized subintervals does not scale because
///        the distribution of the special leaves is highly skewed
///        and most special leaves are in the first few segments
///        whereas later on there are very few special leaves.
///
///        This LoadBalancer gradually increases the number of
///        segments to sieve as long the expected runtime of the
///        sieve distance is smaller than the expected finish time
///        of the algorithm. Near the end the LoadBalancer will
///        gradually decrease the number of segments to sieve in
///        order to prevent that 1 thread will run much longer
///        than all the other threads.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <LoadBalancer.hpp>
#include <primecount-internal.hpp>
#include <json.hpp>
#include <backup.hpp>
#include <S2Status.hpp>
#include <Sieve.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <print.hpp>

#include <stdint.h>
#include <cmath>
#include <string>

using namespace std;

namespace primecount {

int LoadBalancer::get_resume_threads() const
{
  if (is_resume(copy_, "D", x_, y_, z_, k_))
    return calculate_resume_threads(copy_, "D");

  return 0;
}

/// backup resume thread
void LoadBalancer::backup(int thread_id)
{
  double percent = status_.getPercent(low_, sieve_limit_, sum_, sum_approx_);
  double last_backup_seconds = get_time() - backup_time_;

  string tid = "thread" + to_str(thread_id);

  auto& D = json_["D"];
  D["sum"] = to_str(sum_);
  D["percent"] = percent;
  D["seconds"] = get_time() - time_;
  D.erase(tid);

  if (last_backup_seconds > 60)
  {
    backup_time_ = get_time();
    store_backup(json_);
  }
}

/// backup regular thread
void LoadBalancer::backup(int thread_id,
                          int64_t low,
                          int64_t segments,
                          int64_t segment_size)
{
  double percent = status_.getPercent(low_, sieve_limit_, sum_, sum_approx_);
  double last_backup_seconds = get_time() - backup_time_;

  auto& D = json_["D"];
  D["x"] = to_str(x_);
  D["y"] = y_;
  D["z"] = z_;
  D["k"] = k_;
  D["x_star"] = x_star_;
  D["alpha_y"] = get_alpha_y(x_, y_);
  D["alpha_z"] = get_alpha_z(y_, z_);
  D["low"] = low_;
  D["segments"] = segments_;
  D["segment_size"] = segment_size_;
  D["sieve_limit"] = sieve_limit_;
  D["sum"] = to_str(sum_);
  D["percent"] = percent;
  D["seconds"] = get_time() - time_;

  string tid = "thread" + to_str(thread_id);

  if (low <= sieve_limit_)
  {
    D[tid]["low"] = low;
    D[tid]["segments"] = segments;
    D[tid]["segment_size"] = segment_size;
  }
  else
  {
    // finished
    if (D.find(tid) != D.end())
      D.erase(tid);
  }

  if (last_backup_seconds > 60)
  {
    backup_time_ = get_time();
    store_backup(json_);
  }
}

void LoadBalancer::finish_backup()
{
  if (json_.find("D") != json_.end())
    json_.erase("D");

  auto& D = json_["D"];
  D["x"] = to_str(x_);
  D["y"] = y_;
  D["z"] = z_;
  D["k"] = k_;
  D["x_star"] = x_star_;
  D["alpha_y"] = get_alpha_y(x_, y_);
  D["alpha_z"] = get_alpha_z(y_, z_);
  D["sieve_limit"] = sieve_limit_;
  D["sum"] = to_str(sum_);
  D["percent"] = 100.0;
  D["seconds"] = get_time() - time_;

  store_backup(json_);
}

/// resume thread
bool LoadBalancer::resume(int thread_id,
                          int64_t& low,
                          int64_t& segments,
                          int64_t& segment_size) const
{
  if (is_resume(copy_, "D", thread_id, x_, y_, z_, k_))
  {
    string tid = "thread" + to_str(thread_id);
    low = copy_["D"][tid]["low"];
    segments = copy_["D"][tid]["segments"];
    segment_size = copy_["D"][tid]["segment_size"];
    return true;
  }

  return false;
}

// resume result
bool LoadBalancer::resume(maxint_t& sum, double& time) const
{
  if (is_resume(copy_, "D", x_, y_, z_, k_))
  {
    double percent = copy_["D"]["percent"];
    double seconds = copy_["D"]["seconds"];
    print_resume(percent, x_);

    if (!copy_["D"].count("low"))
    {
      sum = to_maxint(copy_["D"]["sum"]);
      time = get_time() - seconds;
      return true;
    }
  }

  return false;
}

LoadBalancer::LoadBalancer(maxint_t x,
                           int64_t y,
                           int64_t z,
                           int64_t k,
                           int64_t sieve_limit,
                           maxint_t sum_approx) :
  x_(x),
  y_(y),
  z_(z),
  k_(k),
  x_star_(get_x_star_gourdon(x, y)),
  low_(0),
  max_low_(0),
  sieve_limit_(sieve_limit),
  segments_(1),
  sum_(0),
  sum_approx_(sum_approx),
  time_(get_time()),
  backup_time_(get_time()),
  status_(x),
  json_(load_backup())
{
  // Read only json copy that is accessed
  // in parallel by multiple threads
  copy_ = json_;

  // start with a tiny segment_size as most
  // special leaves are in the first few segments
  // and we need to ensure that all threads are
  // assigned an equal amount of work
  int64_t sqrt_limit = isqrt(sieve_limit);
  int64_t log = ilog(sqrt_limit);
  log = max(log, 1);
  segment_size_ = sqrt_limit / log;

  int64_t min_size = 1 << 9;
  segment_size_ = max(segment_size_, min_size);
  segment_size_ = Sieve::get_segment_size(segment_size_);

  // try to use a segment size that fits exactly
  // into the CPUs L1 data cache
  int64_t l1_dcache_size = 1 << 15;
  max_size_ = l1_dcache_size * 30;
  max_size_ = max(max_size_, sqrt_limit);
  max_size_ = Sieve::get_segment_size(max_size_);

  if (is_resume(json_, "D", x_, y_, z_, k_))
  {
    double seconds = copy_["D"]["seconds"];
    sum_ = to_maxint(copy_["D"]["sum"]);
    time_ = get_time() - seconds;

    if (json_["D"].count("low"))
    {
      low_ = copy_["D"]["low"];
      segments_ = copy_["D"]["segments"];
      segment_size_ = copy_["D"]["segment_size"];
    }
  }
  else
  {
    if (json_.find("D") != json_.end())
      json_.erase("D");
  }
}

maxint_t LoadBalancer::get_sum() const
{
  return sum_;
}

double LoadBalancer::get_wtime() const
{
  return time_;
}

// Used by resume threads
void LoadBalancer::update_result(int thread_id, uint64_t high, maxint_t sum)
{
  #pragma omp critical (get_work)
  {
    sum_ += sum;
    backup(thread_id);

    if (is_print())
      status_.print(high, sieve_limit_, sum_, sum_approx_);
  }
}

bool LoadBalancer::get_work(int thread_id,
                            int64_t* low,
                            int64_t* segments,
                            int64_t* segment_size,
                            maxint_t sum,
                            Runtime& runtime)
{
  #pragma omp critical (get_work)
  {
    sum_ += sum;

    if (is_print())
    {
      uint64_t dist = *segments * *segment_size;
      uint64_t high = *low + dist;
      status_.print(high, sieve_limit_, sum_, sum_approx_);
    }

    update(low, segments, runtime);

    *low = low_;
    *segments = segments_;
    *segment_size = segment_size_;
    low_ += segments_ * segment_size_;
    low_ = min(low_, sieve_limit_ + 1);

    backup(thread_id, *low, *segments, *segment_size);
  }

  return *low < sieve_limit_;
}

void LoadBalancer::update(int64_t* low,
                          int64_t* segments,
                          Runtime& runtime)
{
  if (*low > max_low_)
  {
    max_low_ = *low;
    segments_ = *segments;

    // We only start increasing the segment_size and segments
    // per thread once the first special leaves have been
    // found. Near the start there is a very large number of
    // leaves and we don't want a single thread to compute
    // them all by himself (which would cause scaling issues).
    if (sum_ == 0)
      return;

    if (segment_size_ < max_size_)
      segment_size_ = min(segment_size_ * 2, max_size_);
    else
      update_segments(runtime);
  }
}

/// Remaining seconds till finished
double LoadBalancer::remaining_secs() const
{
  double percent = status_.getPercent(low_, sieve_limit_, sum_, sum_approx_);
  percent = in_between(10, percent, 100);
  double total_secs = get_time() - time_;
  double secs = total_secs * (100 / percent) - total_secs;
  return secs;
}

/// Increase or decrease the number of segments per thread
/// based on the remaining runtime.
///
void LoadBalancer::update_segments(Runtime& runtime)
{
  // Near the end it is important that threads run only for
  // a short amount of time in order to ensure that all
  // threads finish nearly at the same time. Since the
  // remaining time is just a rough estimation we want to be
  // very conservative so we divide the remaining time by 4.
  double rem_secs = remaining_secs() / 4;

  // If the previous thread runtime is larger than the
  // estimated remaining time the factor that we calculate
  // below will be < 1 and we will reduce the number of
  // segments per thread. Otherwise if the factor > 1 we
  // will increase the number of segments per thread.
  double min_secs = 0.01;
  double divider = max(min_secs, runtime.secs);
  double factor = rem_secs / divider;

  // For small and medium computations the thread runtime
  // should be about 1000x the thread initialization time.
  // However for very large computations we want to further
  // reduce the thread runtimes in order to increase the
  // backup frequency. If the thread runtime is > 6 hours
  // we reduce the thread runtime to about 50x the thread
  // initialization time.
  double init_secs = max(min_secs, runtime.init);
  double init_factor = in_between(50, (3600 * 6) / init_secs, 1000);

  // Reduce the thread runtime if it is much larger than
  // its initialization time. This increases the number of
  // backups without deteriorating performance as we make
  // sure that the thread runtime is still much larger than
  // the thread initialization time.
  if (runtime.secs > min_secs &&
      runtime.secs > runtime.init * init_factor)
  {
    double old = factor;
    double next_runtime = runtime.init * init_factor;
    factor = next_runtime / runtime.secs;
    factor = min(factor, old);
  }

  // Near the end when the remaining time goes close to 0
  // the load balancer tends to reduce the number of
  // segments per thread also close to 0 which is very bad
  // for performance. The condition below fixes this issue
  // and ensures that the thread runtime is always at
  // least 10x the thread initialization time.
  if (runtime.secs > 0 &&
      runtime.secs * factor < runtime.init * 10)
  {
    double next_runtime = runtime.init * 10;
    double current_runtime = runtime.secs;
    factor = next_runtime / current_runtime;
  }

  // Since the distribution of the special leaves is highly
  // skewed (at the beginning) we want to increase the
  // number of segments per thread very slowly. Because if
  // the previously sieved interval contained no special
  // leaves but the next interval contains many special
  // leaves then sieving the next interval might take orders
  // of magnitude more time even if the interval size is
  // identical.
  factor = in_between(0.5, factor, 2.0);
  double next_runtime = runtime.secs * factor;

  if (next_runtime < min_secs)
    segments_ *= 2;
  else
  {
    double new_segments = round(segments_ * factor);
    segments_ = (int64_t) new_segments;
    segments_ = max(segments_, 1);
  }
}

} // namespace
