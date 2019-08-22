///
/// @file  S2_hard_mpi.cpp
/// @brief Calculate the contribution of the hard special leaves using
///        a prime sieve. This is a distributed implementation using
///        MPI (Message Passing Interface) and OpenMP multi-threading.
///
///        Usually the computation of the hard special leaves
///        requires a binary indexed tree a.k.a. Fenwick tree to count
///        the number of unsieved elements in O(log n) time. But it
///        is actually much faster to simply count the number of
///        unsieved elements directly from the sieve array using the
///        POPCNT instruction. Hence this implementation does not use
///        a binary indexed tree.
///
///        This implementation is based on the paper:
///        Tomás Oliveira e Silva, Computing pi(x): the combinatorial
///        method, Revista do DETUA, vol. 4, no. 6, March 2006,
///        pp. 759-768.
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <PiTable.hpp>
#include <FactorTable.hpp>
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
#include <Sieve.hpp>
#include <print.hpp>

#include <stdint.h>
#include <vector>

using namespace std;
using namespace primecount;

namespace {

/// Compute the contribution of the hard special leaves using a
/// segmented sieve. Each thread processes the interval
/// [low, low + segments * segment_size[.
///
/// Note that in the Deleglise-Rivat paper it is suggested to use a
/// segment size of y. In practice however this uses too much memory
/// especially when using multi-threading. Hence we are using a
/// segment size of sqrt(z) as suggested in Xavier Gourdon's paper.
/// In primecount's implementation a segment size of sqrt(z) seems
/// ideal since slightly increasing the segment size decreases
/// performance because of cache misses and slightly decreasing the
/// segment size also decreases performance.
///
template <typename T, typename FactorTable, typename Primes>
T S2_hard_thread(T x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 int64_t low,
                 int64_t segments,
                 int64_t segment_size,
                 FactorTable& factor,
                 PiTable& pi,
                 Primes& primes,
                 Runtime& runtime)
{
  int64_t low1 = max(low, 1);
  int64_t limit = min(low + segments * segment_size, z + 1);
  int64_t max_b = pi[min3(isqrt(x / low1), isqrt(z), y)];
  int64_t pi_sqrty = pi[isqrt(y)];
  T s2_hard = 0;

  if (c > max_b)
    return s2_hard;

  runtime.init_start();
  Sieve sieve(low, segment_size, max_b);
  auto phi = generate_phi(low, max_b, primes, pi);
  runtime.init_stop();

  // Segmented sieve of Eratosthenes
  for (; low < limit; low += segment_size)
  {
    // current segment [low, high[
    int64_t high = min(low + segment_size, limit);
    low1 = max(low, 1);

    // pre-sieve multiples of first c primes
    sieve.pre_sieve(c, low, high);

    int64_t count_low_high = sieve.count((high - 1) - low);
    int64_t b = c + 1;

    // For c + 1 <= b <= pi_sqrty
    // Find all special leaves: n = primes[b] * m
    // which satisfy: mu[m] != 0 && primes[b] < lpf[m] && low <= (x / n) < high
    for (int64_t end = min(pi_sqrty, max_b); b <= end; b++)
    {
      int64_t prime = primes[b];
      T xp = x / prime;
      int64_t xp_div_high = min(fast_div(xp, high), y);
      int64_t min_m = max(xp_div_high, y / prime);
      int64_t max_m = min(fast_div(xp, low1), y);
      int64_t count = 0;
      int64_t start = 0;

      if (prime >= max_m)
        goto next_segment;

      factor.to_index(&min_m);
      factor.to_index(&max_m);

      for (int64_t m = max_m; m > min_m; m--)
      {
        // mu(m) != 0 && prime < lpf(m)
        if (prime < factor.mu_lpf(m))
        {
          int64_t fm = factor.get_number(m);
          int64_t xpm = fast_div64(xp, fm);
          int64_t stop = xpm - low;
          count += sieve.count(start, stop, low, high, count, count_low_high);
          start = stop + 1;
          int64_t phi_xpm = phi[b] + count;
          int64_t mu_m = factor.mu(m);
          s2_hard -= mu_m * phi_xpm;
        }
      }

      phi[b] += count_low_high;
      count_low_high -= sieve.cross_off(b, prime);
    }

    // For pi_sqrty < b <= pi_sqrtz
    // Find all hard special leaves: n = primes[b] * primes[l]
    // which satisfy: low <= (x / n) < high
    for (; b <= max_b; b++)
    {
      int64_t prime = primes[b];
      T xp = x / prime;
      int64_t xp_div_low = min(fast_div(xp, low1), y);
      int64_t xp_div_high = min(fast_div(xp, high), y);
      int64_t l = pi[min(xp_div_low, z / prime)];
      int64_t min_hard = max(xp_div_high, prime);
      int64_t count = 0;
      int64_t start = 0;

      if (prime >= primes[l])
        goto next_segment;

      for (; primes[l] > min_hard; l--)
      {
        int64_t xpq = fast_div64(xp, primes[l]);
        int64_t stop = xpq - low;
        count += sieve.count(start, stop, low, high, count, count_low_high);
        start = stop + 1;
        int64_t phi_xpq = phi[b] + count;
        s2_hard += phi_xpq;
      }

      phi[b] += count_low_high;
      count_low_high -= sieve.cross_off(b, prime);
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

  int64_t max_prime = min(y, z / isqrt(y));
  PiTable pi(max_prime);

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
      // Unsigned integer division is usually slightly
      // faster than signed integer division
      using UT = typename make_unsigned<T>::type;
      s2_hard = S2_hard_thread((UT) x, y, z, c, low, segments, segment_size, factor, pi, primes, runtime);
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
                     int64_t z,
                     T s2_hard_approx)
{
  T s2_hard = 0;
  int slaves = mpi_num_procs() - 1;

  MpiMsg msg;
  MpiLoadBalancer loadBalancer(x, z, s2_hard_approx);
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
      loadBalancer.get_work(&msg);

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
  double time = get_time();

  if (is_mpi_master_proc())
    s2_hard = S2_hard_mpi_master(x, z, s2_hard_approx);
  else
  {
    FactorTable<uint16_t> factor(y, threads);
    int64_t max_prime = min(y, z / isqrt(y));
    auto primes = generate_primes<int32_t>(max_prime);
    S2_hard_slave(x, y, z, c, primes, factor, threads);
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
  double time = get_time();

  if (is_mpi_master_proc())
    s2_hard = S2_hard_mpi_master(x, z, s2_hard_approx);
  else
  {
    // uses less memory
    if (y <= FactorTable<uint16_t>::max())
    {
      FactorTable<uint16_t> factor(y, threads);
      int64_t max_prime = min(y, z / isqrt(y));
      auto primes = generate_primes<uint32_t>(max_prime);
      S2_hard_slave(x, y, z, c, primes, factor, threads);
    }
    else
    {
      FactorTable<uint32_t> factor(y, threads);
      int64_t max_prime = min(y, z / isqrt(y));
      auto primes = generate_primes<int64_t>(max_prime);
      S2_hard_slave(x, y, z, c, primes, factor, threads);
    }
  }

  print("S2_hard", s2_hard, time);
  return s2_hard;
}

#endif

} // namespace
