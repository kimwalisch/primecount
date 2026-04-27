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
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
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
#include <TryLockGuard.hpp>

#include <stdint.h>
#include <atomic>
#include <cmath>

namespace {

// For small PrimePi(x) computations with x ≤ 10^18 the
// best performance is usually achieved using a sieve
// array size that matches the CPU's L1 data cache size
// (per core) or that is slightly larger than the L1 cache
// size but smaller than the L2 cache size (per core).

constexpr int64_t numbers_per_byte = 30;
constexpr int64_t segment_size_alignment = 240;
constexpr int64_t L1_segment_size = L1_CACHE_SIZE * numbers_per_byte;
constexpr int64_t L2_segment_size = L2_CACHE_SIZE * numbers_per_byte;
constexpr int64_t max_packed_segment_size =
    UINT32_MAX - (UINT32_MAX % segment_size_alignment);

} // namespace

namespace primecount {

LoadBalancerS2::LoadBalancerS2(maxint_t x,
                               int64_t y,
                               int64_t sieve_limit,
                               int threads,
                               bool is_print) :
  sieve_limit_(sieve_limit),
  sqrt_limit_(isqrt(sieve_limit)),
  start_time_(get_time()),
  threads_(threads),
  is_print_(is_print),
  status_(x, y, is_print)
{
  int64_t segment_size;
  int64_t segments;

  if (threads == 1 &&
      !is_print)
  {
    segment_size = L1_segment_size;
    segment_size = min(segment_size, sieve_limit);
    segment_size = min(segment_size, max_packed_segment_size);
    segment_size = Sieve::align_segment_size(segment_size);

    // Currently our Sieve.cpp does not rebalance its
    // counters data structure. However, if we process the
    // computation in chunks then the sieve gets recreated
    // for each new chunk which rebalances the counters.
    // Therefore we limit the number of segments here.
    segments = 100;
  }
  else
  {
    // When using multi-threading, it is important to
    // start with a tiny segment size of x^(1/4) as most
    // special leaves are located in the first few segments
    // and as we need to ensure that all threads are
    // assigned an equal amount of work.
    segment_size = isqrt(isqrt(x));
    segment_size = max(segment_size, 1 << 9);
    segment_size = min(segment_size, max_packed_segment_size);
    segment_size = Sieve::align_segment_size(segment_size);
    segments = 1;
  }

  store_packed(segment_size, segments);
}

/// Remaining seconds till finished
double LoadBalancerS2::remaining_secs(int64_t low) const
{
  double percent = status_.getPercent(low, sieve_limit_);
  percent = in_between(10, percent, 100);
  double total_secs = get_time() - start_time_;
  double secs = total_secs * (100 / percent) - total_secs;
  return secs;
}

/// Pack segment_size & segments into a uint64_t,
/// needed for lockfree atomic data access.
///
void LoadBalancerS2::store_packed(uint64_t segment_size,
                                  uint64_t segments)
{
  ASSERT(segments <= UINT32_MAX);
  ASSERT(segment_size <= UINT32_MAX);
  uint64_t packed = segment_size | (segments << 32);
  segment_data_.store(packed, std::memory_order_relaxed);
}

/// Assign new [low, high[ workload to thread.
/// Multiple threads may call get_work() simultaneously, since
/// this function is not protected by a mutex, it must not
/// modify any shared member variables, except the atomic low_
/// max_low_, segment_data_, found_first_leaf_ variables.
///
bool LoadBalancerS2::get_work(ThreadData& thread)
{
  int64_t print_high = 0;
  int64_t max_low = max_low_.load(std::memory_order_relaxed);
  bool found_first_leaf = found_first_leaf_.load(std::memory_order_relaxed);
  uint64_t segment_data = segment_data_.load(std::memory_order_relaxed);
  int64_t segment_size = segment_data & 0xffffffffu;
  int64_t segments = segment_data >> 32;

  if (is_print_ &&
      thread.low >= max_low)
  {
    int64_t dist = thread.segment_size * thread.segments;
    print_high = thread.low + dist;
  }

  if (thread.sum &&
      !found_first_leaf)
  {
    found_first_leaf_.store(true, std::memory_order_relaxed);
    found_first_leaf = true;
  }

  // We only start increasing the segment size and segments
  // per thread once the first special leaf has been found.
  // Most special leaves are located near the start (near y).
  // Hence, we assign tiny work chunks to the threads in
  // this region to avoid load imbalance.
  if (thread.low <= max_low ||
      !found_first_leaf)
  {
    thread.segment_size = segment_size;
    thread.segments = segments;
  }
  else
  {
    max_low_.store(thread.low, std::memory_order_relaxed);
    run_load_balancing(thread, segment_size);
  }

  // The earlier loads are used for heuristic chunk
  // sizing, it is OK if they are slightly outdated. This
  // fetch_add() reserves unique work for this thread.
  int64_t dist = thread.segment_size * thread.segments;
  thread.low = low_.fetch_add(dist, std::memory_order_relaxed);
  thread.sum = 0;
  thread.init_secs = 0;
  thread.secs = 0;

  // The lockfree critical section above should complete
  // as fast as possible. Hence, printing should be done
  // afterwards since it may incur a system call.
  if (print_high)
    print_S2_status(print_high);

  return thread.low < sieve_limit_;
}

void LoadBalancerS2::run_load_balancing(ThreadData& thread,
                                        int64_t segment_size)
{
  // If segment_size < L1_segment_size then slowly increase
  // the segment size until it reaches L1_segment_size.
  if (segment_size < L1_segment_size)
  {
    segment_size += segment_size / 16;
    segment_size = min(segment_size, L1_segment_size);
    segment_size = min(segment_size, max_packed_segment_size);
    segment_size = Sieve::align_segment_size(segment_size);

    store_packed(segment_size, thread.segments);
    thread.segment_size = segment_size;
    return;
  }

  // If segment_size >= L1_segment_size then slowly increase
  // the segment size until it reaches L2_segment_size.
  if (segment_size >= L1_segment_size &&
      segment_size < L2_segment_size &&
      segment_size < sqrt_limit_)
  {
    segment_size += segment_size / 16;
    segment_size = min(segment_size, L2_segment_size);
    segment_size = min(segment_size, max_packed_segment_size);
    segment_size = Sieve::align_segment_size(segment_size);

    store_packed(segment_size, thread.segments);
    thread.segment_size = segment_size;
    return;
  }

  int64_t low = low_.load(std::memory_order_relaxed);
  int64_t segments = get_segments(thread, low);

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
  if (segment_size >= L2_segment_size &&
      segment_size < sqrt_limit_)
  {
    int64_t dist = (segment_size * segments) * threads_;
    int64_t high = min(low + dist, sieve_limit_);

    if (segment_size < isqrt(high))
    {
      segment_size += segment_size / 16;
      dist = (segment_size * segments) * threads_;
      high = min(low + dist, sieve_limit_);
      segment_size = isqrt(high);
      segment_size = Sieve::align_segment_size(segment_size);
    }
  }

  store_packed(segment_size, segments);
  thread.segment_size = segment_size;
  thread.segments = segments;
}

/// Increase or decrease the number of segments per
/// thread based on the remaining runtime.
///
int64_t LoadBalancerS2::get_segments(const ThreadData& thread,
                                     int64_t low) const
{
  // Near the end it is important that threads run only for
  // a short amount of time in order to ensure that all
  // threads finish nearly at the same time. Since the
  // remaining time is just a rough estimation we want to be
  // very conservative so we divide the remaining time by 3.
  double rem_secs = remaining_secs(low) / 3;

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
  int64_t segments = thread.segments;

  if (next_runtime < min_secs)
    segments *= 2;
  else
  {
    double new_segments = std::round(segments * factor);
    segments = (int64_t) new_segments;
    segments = max(segments, 1);
  }

  segments = min(segments, UINT32_MAX);
  return segments;
}

void LoadBalancerS2::print_S2_status(int64_t high)
{
  double time = get_time();

#if __cplusplus >= 201703L
  // Prevent lock contention on many-core systems,
  // print status only every 0.1 seconds.
  if (std::atomic<double>::is_always_lock_free &&
      time <= next_print_time_.load(std::memory_order_relaxed))
    return;
#endif

  // For printing the status it is OK to use a non-blocking
  // userspace lock because printing the status is a non
  // essential operation and hence even if the OS preempts
  // the thread holding the lock it won't cause any deadlocks
  // or performance issues, it will only delay the status
  // output.
  TryLockGuard guard(print_lock_);

  if (guard.owns_lock())
  {
  #if __cplusplus >= 201703L
    // It is theoretically possible that multiple threads
    // enter this critical section within 0.1 seconds.
    // This additional condition prevents it.
    if (std::atomic<double>::is_always_lock_free &&
        time <= next_print_time_.load(std::memory_order_relaxed))
      return;
  #endif

    // The next thread can print again in 0.1 seconds.
    next_print_time_.store(time + 0.1, std::memory_order_relaxed);
    status_.print_S2_hard(high, sieve_limit_);
  }
}

} // namespace
