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
///        This LoadBalancerS2 therefore starts with a tiny segment
///        size of x^(1/4) and one segment per thread. Then the
///        LoadBalancerS2 gradually increases the segment size until
///        it reaches sqrt(limit). Afterwards the LoadBalancerS2
///        gradually increases the number of segments per thread as
///        long as the expected thread runtime is smaller than the
///        expected finish time of the algorithm. Near the end the
///        LoadBalancerS2 gradually decreases the number of segments
///        per thread in order to prevent 1 thread from running much
///        longer than all the other threads.
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

namespace {

// For small PrimePi(x) computations with x ≤ 10^18 the
// best performance is usually achieved using a sieve
// array size that matches the CPU's L1 data cache size
// (per core) or that is slightly larger than the L1 cache
// size but smaller than the L2 cache size (per core).

constexpr int64_t numbers_per_byte = 30;
constexpr int64_t L1_segment_size = L1_CACHE_SIZE * numbers_per_byte;
constexpr int64_t L2_segment_size = L2_CACHE_SIZE * numbers_per_byte;

} // namespace

namespace primecount {

LoadBalancerS2::LoadBalancerS2(maxint_t x,
                               int64_t sieve_limit,
                               maxint_t sum_approx,
                               int threads,
                               bool is_print) :
  sieve_limit_(sieve_limit),
  sqrt_limit_(isqrt(sieve_limit)),
  sum_approx_(sum_approx),
  time_(get_time()),
  threads_(threads),
  is_print_(is_print),
  status_(x)
{
  lock_.init(threads);

  if (threads == 1 &&
      !is_print)
  {
    segment_size_ = L1_segment_size;
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
  std::string status;
  bool has_work;

  {
    LockGuard lockGuard(lock_);
    sum_ += thread.sum;

    if (is_print_)
    {
      uint64_t dist = thread.segment_size * thread.segments;
      uint64_t high = thread.low + dist;
      status = status_.getStatus(high, sieve_limit_, sum_, sum_approx_);
    }

    update_load_balancing(thread);

    thread.low = low_;
    thread.segments = segments_;
    thread.segment_size = segment_size_;
    thread.sum = 0;
    thread.secs = 0;
    thread.init_secs = 0;

    has_work = thread.low < sieve_limit_;
    low_ += segment_size_ * segments_;
  }

  // Printing to the terminal incurs a system call
  // and may hence be slow. Therefore, we do it
  // after having released the mutex.
  if (!status.empty())
    print_status(status);

  return has_work;
}

void LoadBalancerS2::update_load_balancing(const ThreadData& thread)
{
  if (thread.low > max_low_)
  {
    max_low_ = thread.low;
    segments_ = thread.segments;

    // We only start increasing the segment size and segments per
    // thread once the first special leaves have been found. Most
    // special leaves are located near the start (around y).
    // Hence, we assign tiny work chunks to the threads in this
    // region to avoid load imbalance.
    if (sum_ == 0)
      return;

    // If segment_size < L1_segment_size then slowly increase the
    // segment size until it reaches L1_segment_size.
    if (segment_size_ < L1_segment_size)
    {
      segment_size_ += segment_size_ / 16;
      segment_size_ = min(segment_size_, L1_segment_size);
      segment_size_ = Sieve::align_segment_size(segment_size_);
      return;
    }

    // If segment_size >= L1_segment_size then slowly increase the
    // segment size until it reaches L2_segment_size.
    if (segment_size_ >= L1_segment_size &&
        segment_size_ < L2_segment_size &&
        segment_size_ < sqrt_limit_)
    {
      segment_size_ += segment_size_ / 16;
      segment_size_ = min(segment_size_, L2_segment_size);
      segment_size_ = Sieve::align_segment_size(segment_size_);
      return;
    }

    update_number_of_segments(thread);

    // The hard special leaves algorithm is basically a modified
    // segmented sieve of Eratosthenes. Using the segmented sieve of
    // Eratosthenes it is of utmost importance that sieve array fits
    // into the CPU's cache, otherwise performance will deteriorate
    // significantly and the algorithm will scale poorly.
    //
    // Deleglise-Rivat orignially suggested using a segment size of
    // O(y). Xavier Gourdon realized this segment size was much too
    // large for new record PrimePi(x) computations and hence
    // suggested using a smaller segment size of O(sqrt(x/y)) which
    // is the same as O(sqrt(sieve_limit)). However, for new record
    // PrimePi(x) computations a segment size O(sqrt(sieve_limit))
    // is still too large. Hence, I use an even smaller segment size
    // of O(sqrt(high)) in primecount.
    if (segment_size_ >= L2_segment_size &&
        segment_size_ < sqrt_limit_)
    {
      int64_t dist = (segment_size_ * segments_) * threads_;
      int64_t high = min(low_ + dist, sieve_limit_);

      if (segment_size_ < isqrt(high))
      {
        segment_size_ += segment_size_ / 16;
        dist = (segment_size_ * segments_) * threads_;
        high = min(low_ + dist, sieve_limit_);
        segment_size_ = isqrt(high);
        segment_size_ = Sieve::align_segment_size(segment_size_);
      }
    }
  }
}

/// Increase or decrease the number of segments per
/// thread based on the remaining runtime.
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
