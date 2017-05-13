///
/// @file  pi_lmo_parallel.cpp
/// @brief Parallel implementation of the Lagarias-Miller-Odlyzko
///        prime counting algorithm using OpenMP. This implementation
///        uses load balancing and counts the number of unsieved
///        elements using POPCNT without using any special counting
///        tree data structure.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <BitSieve.hpp>
#include <generate.hpp>
#include <generate_phi.hpp>
#include <min.hpp>
#include <imath.hpp>
#include <PhiTiny.hpp>
#include <PiTable.hpp>
#include <S1.hpp>
#include <S2Status.hpp>
#include <S2LoadBalancer.hpp>
#include <Wheel.hpp>

#include <stdint.h>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;

namespace {

/// Cross-off the multiples of prime in the sieve array
int64_t cross_off(BitSieve& sieve,
                  int64_t low,
                  int64_t high,
                  int64_t prime,
                  WheelItem& w)
{
  int64_t unset = 0;
  int64_t m = w.next_multiple;
  int64_t wheel_index = w.wheel_index;

  for (; m < high; m += prime * Wheel::next_multiple_factor(&wheel_index))
  {
    // +1 if m is unset the first time
    unset += sieve[m - low];
    sieve.unset(m - low);
  }

  w.set(m, wheel_index);
  return unset;
}

/// Compute the S2 contribution for the interval
/// [low_process, low_process + segments * segment_size[.
/// The missing special leaf contributions for the interval
/// [1, low_process[ are later reconstructed and added in
/// the calling (parent) S2 function.
///
int64_t S2_thread(int64_t x,
                  int64_t y,
                  int64_t c,
                  int64_t segment_size,
                  int64_t segments_per_thread,
                  int64_t low,
                  int64_t limit,
                  PiTable& pi,
                  vector<int32_t>& primes,
                  vector<int32_t>& lpf,
                  vector<int32_t>& mu)
{
  limit = min(low + segment_size * segments_per_thread, limit);
  int64_t size = pi[min(isqrt(x / low), y)] + 1;
  int64_t pi_sqrty = pi[isqrt(y)];
  int64_t pi_y = pi[y];
  int64_t S2_thread = 0;

  if (c >= size - 1)
    return 0;

  BitSieve sieve(segment_size);
  Wheel wheel(primes, size, low);
  auto phi = generate_phi(low - 1, size - 1, primes, pi);

  // segmented sieve of Eratosthenes
  for (; low < limit; low += segment_size)
  {
    // current segment = [low, high[
    int64_t high = min(low + segment_size, limit);
    int64_t b = c + 1;

    // pre-sieve the multiples of the first c primes
    sieve.pre_sieve(c, low);

    int64_t count_low_high = sieve.count((high - 1) - low);

    // For c + 1 <= b < pi_sqrty
    // Find all special leaves: n = primes[b] * m
    // which satisfy:  mu[m] != 0 && primes[b] < lpf[m], low <= (x / n) < high
    for (; b < min(pi_sqrty, size); b++)
    {
      int64_t prime = primes[b];
      int64_t min_m = max(x / (prime * high), y / prime);
      int64_t max_m = min(x / (prime * low), y);
      int64_t count = 0;
      int64_t i = 0;

      if (prime >= max_m)
        goto next_segment;

      for (int64_t m = max_m; m > min_m; m--)
      {
        if (mu[m] != 0 && prime < lpf[m])
        {
          int64_t xn = x / (prime * m);
          int64_t stop = xn - low;
          count += sieve.count(i, stop);
          i = stop + 1;
          int64_t phi_xn = phi[b] + count;
          S2_thread -= mu[m] * phi_xn;
        }
      }

      phi[b] += count_low_high;
      count_low_high -= cross_off(sieve, low, high, prime, wheel[b]);
    }

    // For pi_sqrty <= b < pi_y
    // Find all special leaves: n = primes[b] * prime2
    // which satisfy: low <= (x / n) < high
    for (; b < min(pi_y, size); b++)
    {
      int64_t prime = primes[b];
      int64_t l = pi[min(x / (prime * low), y)];
      int64_t min_m = max(x / (prime * high), y / prime, prime);
      int64_t count = 0;
      int64_t i = 0;

      if (prime >= primes[l])
        goto next_segment;

      for (; primes[l] > min_m; l--)
      {
        int64_t xn = x / (prime * primes[l]);
        int64_t stop = xn - low;
        count += sieve.count(i, stop);
        i = stop + 1;
        int64_t phi_xn = phi[b] + count;
        S2_thread += phi_xn;
      }

      phi[b] += count_low_high;
      count_low_high -= cross_off(sieve, low, high, prime, wheel[b]);
    }

    next_segment:;
  }

  return S2_thread;
}

class NewLoadBalancer
{
public:
  NewLoadBalancer(maxint_t x,
                  int64_t y,
                  int64_t limit,
                  maxint_t s2_approx,
                  int threads)
    : status_(x),
      loadBalancer_(x, y, limit, threads),
      low_(1),
      limit_(x / y + 1),
      segments_(1),
      max_size_(next_power_of_2(max(isqrt(x / y), 1024))),
      S2_total_(0),
      time_(get_wtime()),
      s2_approx_(s2_approx),
      finished_(false)
  {
    segment_size_ = loadBalancer_.get_min_segment_size();
  }

  bool is_increase(double percent, double seconds, double elapsed_time)
  {
    double min_seconds = 0.01;

    if (seconds < min_seconds)
      return true;

    // avoid division by 0
    percent = in_between(1, percent, 99);

    // calculate remaining time till finished
    double remaining_time = elapsed_time * (100 / percent) - elapsed_time;
    double max_seconds = remaining_time / 6;
    double is_increase = max(min_seconds, max_seconds);

    return seconds < is_increase;
  }

  void update(int64_t* low, int64_t* segments, int64_t* segment_size, maxint_t S2, double seconds)
  {
    double elapsed_time = get_wtime() - time_;

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

        if (is_increase(percent, seconds, elapsed_time))
        {
          if (segment_size_ < max_size_)
            segment_size_ *= 2;
          else
            segments_ += (segments_ / 2) + 1;
        }
        else
          segments_ -= segments_ / 3;
      }
    }

    if (is_print())
      status_.print(S2_total_, s2_approx_);
  }
  
  bool finished()
  {
    return finished_;
  }

  maxint_t getResult()
  {
    return S2_total_;
  }

private:
  S2Status status_;
  S2LoadBalancer loadBalancer_;
  int64_t low_;
  int64_t limit_;
  int64_t segment_size_;
  int64_t segments_;
  int64_t max_size_;
  maxint_t S2_total_;
  double time_;
  maxint_t s2_approx_;
  bool finished_;
};

/// Calculate the contribution of the special leaves.
/// This is a parallel implementation with advanced load balancing.
/// As most special leaves tend to be in the first segments we
/// start off with a small segment size and few segments
/// per thread, after each iteration we dynamically increase
/// the segment size and the segments per thread.
///
int64_t S2(int64_t x,
           int64_t y,
           int64_t c,
           int64_t s2_approx,
           vector<int32_t>& primes,
           vector<int32_t>& lpf,
           vector<int32_t>& mu,
           int threads)
{
  print("");
  print("=== S2(x, y) ===");
  print("Computation of the special leaves");

  int64_t limit = x / y + 1;
  double time = get_wtime();
  threads = ideal_num_threads(threads, limit);
  NewLoadBalancer loadBalancer(x, y, limit, s2_approx, threads);
  PiTable pi(y);

  #pragma omp parallel for num_threads(threads)
  for (int i = 0; i < threads; i++)
  {
    int64_t low;
    int64_t segments;
    int64_t segment_size;
    int64_t S2 = 0;
    double seconds = 0;

    while (true)
    {
      loadBalancer.update(&low, &segments, &segment_size, S2, seconds);

      if (loadBalancer.finished())
        break;

      seconds = get_wtime();
      S2 = S2_thread(x, y, c, segment_size, segments, low, limit, pi, primes, lpf, mu);
      seconds = get_wtime() - seconds;
    }
  }

  int64_t S2_total = loadBalancer.getResult();
  print("S2", S2_total, time);

  return S2_total;
}

} // namespace

namespace primecount {

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3) / log x)
/// Memory usage: O(x^(1/3) * (log x)^2)
///
int64_t pi_lmo_parallel(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  double alpha = get_alpha_lmo(x);
  int64_t x13 = iroot<3>(x);
  int64_t y = (int64_t) (x13 * alpha);
  int64_t z = x / y;
  int64_t c = PhiTiny::get_c(y);

  print("");
  print("=== pi_lmo_parallel(x) ===");
  print("pi(x) = S1 + S2 + pi(y) - 1 - P2");
  print(x, y, z, c, alpha, threads);

  int64_t p2 = P2(x, y, threads);
  auto primes = generate_primes<int32_t>(y);
  auto lpf = generate_lpf(y);
  auto mu = generate_moebius(y);

  int64_t pi_y = primes.size() - 1;
  int64_t s1 = S1(x, y, c, threads);
  int64_t s2_approx = S2_approx(x, pi_y, p2, s1);
  int64_t s2 = S2(x, y, c, s2_approx, primes, lpf, mu, threads);
  int64_t phi = s1 + s2;
  int64_t sum = phi + pi_y - 1 - p2;

  return sum;
}

} // namespace
