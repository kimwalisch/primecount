///
/// @file  LoadBalancerS2.cpp
/// @brief The LoadBalancerS2 assigns work to the individual
///        threads in the computation of the special leaves in the
///        Lagarias-Miller-Odlyzko, Deleglise-Rivat and Gourdon
///        prime counting algorithms. This load balancer is used by
///        the S2_hard(x, y) and D(x, y) functions.
///
///        Simply parallelizing the computation of the special leaves
///        in the Lagarias-Miller-Odlyzko algorithm by subdividing
///        the sieve interval by the number of threads into equally
///        sized subintervals does not scale because the distribution
///        of the special leaves is highly skewed and most special
///        leaves are in the first few segments whereas later on
///        there are very few special leaves.
///
///        This LoadBalancerS2 gradually increases the number of
///        segments to sieve as long as the expected runtime of the
///        sieve distance is smaller than the expected finish time
///        of the algorithm. Near the end the LoadBalancerS2 will
///        gradually decrease the number of segments to sieve in
///        order to prevent that 1 thread will run much longer than
///        all the other threads.
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <LoadBalancerS2.hpp>
#include <primecount-config.hpp>
#include <primecount-internal.hpp>
#include <StatusS2.hpp>
#include <Sieve.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <min.hpp>

#include <stdint.h>

namespace primecount {

LoadBalancerS2::LoadBalancerS2(maxint_t x,
                               int64_t sieve_limit,
                               maxint_t sum_approx,
                               int threads,
                               bool is_print) :
  sieve_limit_(sieve_limit),
  sum_approx_(sum_approx),
  time_(get_time()),
  is_print_(is_print),
  status_(x)
{
  lock_.init(threads);

  // Using a single thread, the best performance is
  // usually achieved using a sieve array size that matches
  // your CPU's L1 data cache size (per core) or that is
  // slightly larger than your L1 cache size but smaller
  // than your L2 cache size (per core).
  int64_t numbers_per_byte = 30;
  int64_t cache_bytes = L1D_CACHE_SIZE * 2;
  cache_segment_size_ = cache_bytes * numbers_per_byte;
  cache_segment_size_ = Sieve::align_segment_size(cache_segment_size_);

  if (threads == 1 &&
      !is_print)
  {
    segment_size_ = cache_segment_size_;
    segment_size_ = min(segment_size_, sieve_limit);
    segment_size_ = Sieve::align_segment_size(segment_size_);

    // Currently our Sieve.cpp does not rebalance its
    // counters data structure. However, if we process the
    // computation in chunks then the sieve gets recreated
    // for each new chunk which rebalances the counters.
    // Therefore we limit the number of segments here.
    segments_ = 100;
  }
  else
  {
    // When using multi-threading, it is important to
    // start with a tiny segment size of x^(1/4) as most
    // special leaves are located in the first few segments
    // and as we need to ensure that all threads are
    // assigned an equal amount of work.
    segment_size_ = isqrt(isqrt(x));
    segment_size_ = max(segment_size_, 1 << 9);
    segment_size_ = Sieve::align_segment_size(segment_size_);
    segments_ = 1;
  }
}

maxint_t LoadBalancerS2::get_sum() const
{
  return sum_;
}

bool LoadBalancerS2::get_work(ThreadData& thread)
{
  LockGuard lockGuard(lock_);
  sum_ += thread.sum;

  if (is_print_)
  {
    uint64_t dist = thread.segment_size * thread.segments;
    uint64_t high = thread.low + dist;
    status_.print(high, sieve_limit_, sum_, sum_approx_);
  }

  update_load_balancing(thread);

  thread.low = low_;
  thread.segments = segments_;
  thread.segment_size = segment_size_;
  thread.sum = 0;
  thread.secs = 0;
  thread.init_secs = 0;

  low_ += segment_size_ * segments_;
  bool is_work = thread.low < sieve_limit_;

  return is_work;
}

void LoadBalancerS2::update_load_balancing(const ThreadData& thread)
{
  if (thread.low > max_low_)
  {
    max_low_ = thread.low;
    segments_ = thread.segments;

    // We only start increasing the segment_size and segments
    // per thread once the first special leaves have been
    // found. Near the start there is a very large number of
    // leaves and we don't want a single thread to compute
    // them all by himself (which would cause scaling issues).
    if (sum_ == 0)
      return;

    if (segment_size_ >= cache_segment_size_)
      update_number_of_segments(thread);

    // The segmented sieve of Eratosthenes traditionally
    // requires using a segment size of O(sqrt(x)), using a
    // smaller segment size deteriorates the runtime complexity.
    // However, it is possible to use a smaller segment size
    // of O(sqrt(high)) provided that it is dynamically
    // increased after each sieved sieved segment. Using this
    // smaller segment size of O(sqrt(high)) uses less memory
    // and is hence more cache efficient.
    int64_t high = low_ + segment_size_ * segments_;
    high = min(high, sieve_limit_);
    int64_t new_segment_size = isqrt(high);
    new_segment_size = max(cache_segment_size_, new_segment_size);

    // Slowly increase the segment size until it reaches
    // sqrt(high). Most special leaves are located around y,
    // hence we need to be careful to not assign too much work
    // (i.e. use too large segment size) to a single thread
    // in this region.
    if (segment_size_ < new_segment_size)
    {
      segment_size_ += segment_size_ / 16;
      segment_size_ = min(segment_size_, new_segment_size);
      segment_size_ = Sieve::align_segment_size(segment_size_);
    }
  }
}

/// Increase or decrease the number of segments per thread
/// based on the remaining runtime.
///
void LoadBalancerS2::update_number_of_segments(const ThreadData& thread)
{
  // Near the end it is important that threads run only for
  // a short amount of time in order to ensure that all
  // threads finish nearly at the same time. Since the
  // remaining time is just a rough estimation we want to be
  // very conservative so we divide the remaining time by 3.
  double rem_secs = remaining_secs() / 3;

  // If the previous thread runtime is larger than the
  // estimated remaining time the factor that we calculate
  // below will be < 1 and we will reduce the number of
  // segments per thread. Otherwise if the factor > 1 we
  // will increase the number of segments per thread.
  double min_secs = 0.001;
  double divider = max(min_secs, thread.secs);
  double factor = rem_secs / divider;

  // For small and medium computations the thread runtime
  // should be about 5000x the thread initialization time.
  // However for very large computations we want to further
  // reduce the thread runtimes in order to increase the
  // backup frequency. If the thread runtime is > 6 hours
  // we reduce the thread runtime to about 200x the thread
  // initialization time.
  double init_secs = max(min_secs, thread.init_secs);
  double init_factor = in_between(200, (3600 * 6) / init_secs, 5000);

  // Reduce the thread runtime if it is much larger than
  // its initialization time. This increases the number of
  // backups without deteriorating performance as we make
  // sure that the thread runtime is still much larger than
  // the thread initialization time.
  if (thread.secs > min_secs &&
      thread.secs > thread.init_secs * init_factor)
  {
    double old = factor;
    double next_runtime = thread.init_secs * init_factor;
    factor = next_runtime / thread.secs;
    factor = min(factor, old);
  }

  // Near the end when the remaining time goes close to 0
  // the load balancer tends to reduce the number of
  // segments per thread also close to 0 which is very bad
  // for performance. The condition below fixes this issue
  // and ensures that the thread runtime is always at
  // least 20x the thread initialization time.
  if (thread.secs > 0 &&
      thread.secs * factor < thread.init_secs * 20)
  {
    double next_runtime = thread.init_secs * 20;
    double current_runtime = thread.secs;
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
  double next_runtime = thread.secs * factor;

  if (next_runtime < min_secs)
    segments_ *= 2;
  else
  {
    double new_segments = std::round(segments_ * factor);
    segments_ = (int64_t) new_segments;
    segments_ = max(segments_, 1);
  }
}

/// Remaining seconds till finished
double LoadBalancerS2::remaining_secs() const
{
  double percent = status_.getPercent(low_, sieve_limit_, sum_, sum_approx_);
  percent = in_between(10, percent, 100);
  double total_secs = get_time() - time_;
  double secs = total_secs * (100 / percent) - total_secs;
  return secs;
}

} // namespace
