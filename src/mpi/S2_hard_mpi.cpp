///
/// @file  S2_hard_mpi.cpp
/// @brief Calculate the contribution of the hard special leaves
///        (Deleglise-Rivat algorithm) using a segmented sieve.
///        This is a parallel implementation which uses compression
///        (PiTable & FactorTable) to reduce the memory usage by
///        about 10x.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <PiTable.hpp>
#include <FactorTable.hpp>
#include <BinaryIndexedTree.hpp>
#include <BitSieve.hpp>
#include <fast_div.hpp>
#include <generate.hpp>
#include <generate_phi.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <mpi_reduce_sum.hpp>
#include <MpiLoadBalancer.hpp>
#include <MpiMsg.hpp>
#include <imath.hpp>
#include <S2Status.hpp>
#include <S2.hpp>
#include <Wheel.hpp>

#include <stdint.h>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;

namespace {

/// Cross-off the multiples of prime in the sieve array.
/// @return  Count of crossed-off multiples.
///
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

/// Cross-off the multiples of prime in the sieve array.
/// For each element that is unmarked the first time update
/// the binary indexed tree data structure.
///
void cross_off(BitSieve& sieve,
               int64_t low,
               int64_t high,
               int64_t prime,
               WheelItem& w,
               BinaryIndexedTree& tree)
{
  int64_t m = w.next_multiple;
  int64_t wheel_index = w.wheel_index;

  for (; m < high; m += prime * Wheel::next_multiple_factor(&wheel_index))
  {
    if (sieve[m - low])
    {
      sieve.unset(m - low);
      tree.update(m - low);
    }
  }

  w.set(m, wheel_index);
}

/// Returns true if the interval [low, high]
/// contains few hard special leaves
///
bool few_leaves(int64_t low,
                int64_t high,
                int64_t y,
                double alpha)
{
  return (high < y || low > y * alpha);
}

/// Compute the S2 contribution of the hard special leaves
/// using a sieve. Each thread processes the interval
/// [low, low + segments * segment_size[
///
template <typename T, typename FactorTable, typename Primes>
T S2_hard_thread(T x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 int64_t low,
                 int64_t limit,
                 int64_t segments,
                 int64_t segment_size,
                 double alpha,
                 FactorTable& factor,
                 PiTable& pi,
                 Primes& primes,
                 Runtime& runtime)
{
  limit = min(low + segments * segment_size, limit);
  int64_t max_b = pi[min(isqrt(x / low), isqrt(z), y)];
  int64_t pi_sqrty = pi[isqrt(y)];
  T s2_hard = 0;

  if (c > max_b)
    return s2_hard;

  runtime.init_start();
  BinaryIndexedTree tree;
  BitSieve sieve(segment_size);
  Wheel wheel(primes, max_b + 1, low);
  auto phi = generate_phi(low - 1, max_b, primes, pi);
  runtime.init_stop();

  // Segmented sieve of Eratosthenes
  for (; low < limit; low += segment_size)
  {
    // Current segment = interval [low, high[
    int64_t high = min(low + segment_size, limit);
    int64_t b = c + 1;

    // pre-sieve the multiples of the first c primes
    sieve.pre_sieve(c, low);

    // If there are relatively few hard special leaves per segment
    // we count the number of unsieved elements directly from the
    // sieve array using the POPCNT instruction.
    if (few_leaves(low, high, y, alpha))
    {
      int64_t count_low_high = sieve.count((high - 1) - low);

      // For c + 1 <= b <= pi_sqrty
      // Find all special leaves: n = primes[b] * m
      // which satisfy: mu[m] != 0 && primes[b] < lpf[m] && low <= (x / n) < high
      for (int64_t end = min(pi_sqrty, max_b); b <= end; b++)
      {
        int64_t prime = primes[b];
        T x2 = x / prime;
        int64_t x2_div_high = min(fast_div(x2, high), y);
        int64_t min_m = max(x2_div_high, y / prime);
        int64_t max_m = min(fast_div(x2, low), y);
        int64_t count = 0;
        int64_t start = 0;

        if (prime >= max_m)
          goto next_segment;

        factor.to_index(&min_m);
        factor.to_index(&max_m);

        for (int64_t m = max_m; m > min_m; m--)
        {
          if (prime < factor.lpf(m))
          {
            int64_t fm = factor.get_number(m);
            int64_t xn = (int64_t) fast_div(x2, fm);
            int64_t stop = xn - low;
            count += sieve.count(start, stop, low, high, count, count_low_high);
            start = stop + 1;
            int64_t phi_xn = phi[b] + count;
            int64_t mu_m = factor.mu(m);
            s2_hard -= mu_m * phi_xn;
          }
        }

        phi[b] += count_low_high;
        count_low_high -= cross_off(sieve, low, high, prime, wheel[b]);
      }

      // For pi_sqrty < b <= pi_sqrtz
      // Find all hard special leaves: n = primes[b] * primes[l]
      // which satisfy: low <= (x / n) < high
      for (; b <= max_b; b++)
      {
        int64_t prime = primes[b];
        T x2 = x / prime;
        int64_t x2_div_low = min(fast_div(x2, low), y);
        int64_t x2_div_high = min(fast_div(x2, high), y);
        int64_t l = pi[min(x2_div_low, z / prime)];
        int64_t min_hard = max(x2_div_high, prime);
        int64_t count = 0;
        int64_t start = 0;

        if (prime >= primes[l])
          goto next_segment;

        for (; primes[l] > min_hard; l--)
        {
          int64_t xn = (int64_t) fast_div(x2, primes[l]);
          int64_t stop = xn - low;
          count += sieve.count(start, stop, low, high, count, count_low_high);
          start = stop + 1;
          int64_t phi_xn = phi[b] + count;
          s2_hard += phi_xn;
        }

        phi[b] += count_low_high;
        count_low_high -= cross_off(sieve, low, high, prime, wheel[b]);
      }
    }
    else
    {
      // Calculate the contribution of the hard special leaves using
      // Tomás Oliveira's binary indexed tree data structure
      // for counting the number of unsieved elements. This algorithm
      // runs fastest if there are many special leaves per segment.

      // initialize binary indexed tree from sieve
      tree.init(sieve);

      // For c + 1 <= b <= pi_sqrty
      // Find all special leaves: n = primes[b] * m
      // which satisfy: mu[m] != 0 && primes[b] < lpf[m] && low <= (x / n) < high
      for (int64_t end = min(pi_sqrty, max_b); b <= end; b++)
      {
        int64_t prime = primes[b];
        T x2 = x / prime;
        int64_t x2_div_high = min(fast_div(x2, high), y);
        int64_t min_m = max(x2_div_high, y / prime);
        int64_t max_m = min(fast_div(x2, low), y);

        if (prime >= max_m)
          goto next_segment;

        factor.to_index(&min_m);
        factor.to_index(&max_m);

        for (int64_t m = max_m; m > min_m; m--)
        {
          if (prime < factor.lpf(m))
          {
            int64_t fm = factor.get_number(m);
            int64_t xn = (int64_t) fast_div(x2, fm);
            int64_t count = tree.count(low, xn);
            int64_t phi_xn = phi[b] + count;
            int64_t mu_m = factor.mu(m);
            s2_hard -= mu_m * phi_xn;
          }
        }

        phi[b] += tree.count(low, high - 1);
        cross_off(sieve, low, high, prime, wheel[b], tree);
      }

      // For pi_sqrty < b <= pi_sqrtz
      // Find all hard special leaves: n = primes[b] * primes[l]
      // which satisfy: low <= (x / n) < high
      for (; b <= max_b; b++)
      {
        int64_t prime = primes[b];
        T x2 = x / prime;
        int64_t x2_div_low = min(fast_div(x2, low), y);
        int64_t x2_div_high = min(fast_div(x2, high), y);
        int64_t l = pi[min(x2_div_low, z / prime)];
        int64_t min_hard = max(x2_div_high, prime);

        if (prime >= primes[l])
          goto next_segment;

        for (; primes[l] > min_hard; l--)
        {
          int64_t xn = (int64_t) fast_div(x2, primes[l]);
          int64_t count = tree.count(low, xn);
          int64_t phi_xn = phi[b] + count;
          s2_hard += phi_xn;
        }

        phi[b] += tree.count(low, high - 1);
        cross_off(sieve, low, high, prime, wheel[b], tree);
      }
    }

    next_segment:;
  }

  return s2_hard;
}

/// S2_hard MPI slave process.
/// Asks MPI master process for new work and reports
/// partial results to MPI master process.
///
template <typename T, typename FactorTable, typename Primes>
void S2_hard_slave(T x,
                   int64_t y,
                   int64_t z,
                   int64_t c,
                   Primes& primes,
                   FactorTable& factor,
                   int threads)
{
  threads = ideal_num_threads(threads, z);

  int64_t limit = z + 1;
  int64_t max_prime = z / isqrt(y);
  PiTable pi(max_prime);
  double alpha = get_alpha(x, y);

  MpiMsg msg;
  int master_proc_id = mpi_master_proc_id();
  int proc_id = mpi_proc_id();

  #pragma omp parallel for num_threads(threads)
  for (int i = 0; i < threads; i++)
  {
    T s2_hard = 0;
    int64_t low = 0;
    int64_t segments = 0;
    int64_t segment_size = 0;
    Runtime runtime;

    while (true)
    {
      #pragma omp critical (mpi_sync)
      {
        // send result to master process
        msg.set(proc_id, i, low, segments, segment_size, s2_hard, runtime.init, runtime.secs);
        msg.send(master_proc_id);

        // receive new work todo
        msg.recv(proc_id);
        low = msg.low();
        segments = msg.segments();
        segment_size = msg.segment_size();
      }

      if (low > z)
        break;

      runtime.start();
      s2_hard = S2_hard_thread(x, y, z, c, low, limit, segments, segment_size, alpha, factor, pi, primes, runtime);
      runtime.stop();
    }
  }

  msg.set_finished();
  msg.send(master_proc_id);
}

/// S2_hard MPI master process.
/// Assigns work to the MPI slave processes.
///
template <typename T>
T S2_hard_mpi_master(T x,
                     int64_t y,
                     int64_t z,
                     T s2_hard_approx)
{
  T s2_hard = 0;
  int slaves = mpi_num_procs() - 1;

  MpiMsg msg;
  MpiLoadBalancer loadBalancer(x, y, z, s2_hard_approx);
  S2Status status(x);

  while (slaves > 0)
  {
    // wait for results from slave process
    msg.recv_any();

    if (msg.finished())
      slaves--;
    else
    {
      s2_hard += msg.s2_hard<T>();

      // update msg with new work
      loadBalancer.get_work(&msg, s2_hard);

      // send new work to slave process
      msg.send(msg.proc_id());

      if (is_print())
        status.print(s2_hard, s2_hard_approx);
    }
  }

  return s2_hard;
}

} // namespace

namespace primecount {

int64_t S2_hard_mpi(int64_t x,
                    int64_t y,
                    int64_t z,
                    int64_t c,
                    int64_t s2_hard_approx,
                    int threads)
{
  print("");
  print("=== S2_hard_mpi(x, y) ===");
  print("Computation of the hard special leaves");
  print(x, y, c, threads);

  int64_t s2_hard = 0;
  double time = get_wtime();

  if (is_mpi_master_proc())
    s2_hard = S2_hard_mpi_master(x, y, z, s2_hard_approx);
  else
  {
    FactorTable<uint16_t> factor(y, threads);
    int64_t max_prime = z / isqrt(y);
    auto primes = generate_primes<int32_t>(max_prime);

    S2_hard_slave((intfast64_t) x, y, z, c, primes, factor, threads);
  }

  print("S2_hard", s2_hard, time);
  return s2_hard;
}

#ifdef HAVE_INT128_T

int128_t S2_hard_mpi(int128_t x,
                     int64_t y,
                     int64_t z,
                     int64_t c,
                     int128_t s2_hard_approx,
                     int threads)
{
  print("");
  print("=== S2_hard_mpi(x, y) ===");
  print("Computation of the hard special leaves");
  print(x, y, c, threads);

  int128_t s2_hard = 0;
  double time = get_wtime();

  if (is_mpi_master_proc())
    s2_hard = S2_hard_mpi_master(x, y, z, s2_hard_approx);
  else
  {
    // uses less memory
    if (y <= FactorTable<uint16_t>::max())
    {
      FactorTable<uint16_t> factor(y, threads);
      int64_t max_prime = z / isqrt(y);
      auto primes = generate_primes<uint32_t>(max_prime);

      S2_hard_slave((intfast128_t) x, y, z, c, primes, factor, threads);
    }
    else
    {
      FactorTable<uint32_t> factor(y, threads);
      int64_t max_prime = z / isqrt(y);
      auto primes = generate_primes<int64_t>(max_prime);

      S2_hard_slave((intfast128_t) x, y, z, c, primes, factor, threads);
    }
  }

  print("S2_hard", s2_hard, time);
  return s2_hard;
}

#endif

} // namespace
