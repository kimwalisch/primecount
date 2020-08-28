///
/// @file  D_mpi.cpp
/// @brief Implementation of the D formula (from Xavier Gourdon's
///        algorithm) that has been distributed using MPI (Message
///        Passing Interface) and multi-threaded using OpenMP.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <DFactorTable.hpp>
#include <PiTable.hpp>
#include <Sieve.hpp>
#include <LoadBalancer.hpp>
#include <fast_div.hpp>
#include <generate.hpp>
#include <generate_phi.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <print.hpp>
#include <mpi_reduce_sum.hpp>
#include <MpiLoadBalancer.hpp>

#include <stdint.h>
#include <vector>

using namespace std;
using namespace primecount;

namespace {

/// Compute the contribution of the hard special leaves using a
/// segmented sieve. Each thread processes the interval
/// [low, low + segments * segment_size[.
///
template <typename T,
          typename Primes,
          typename DFactorTable,
          typename PhiCache>
T D_thread(T x,
           int64_t x_star,
           int64_t xz,
           int64_t y,
           int64_t z,
           int64_t k,
           const Primes& primes,
           const PiTable& pi,
           const DFactorTable& factor,
           PhiCache& phiCache,
           ThreadSettings& thread)
{
  T sum = 0;

  int64_t low = thread.low;
  int64_t segments = thread.segments;
  int64_t segment_size = thread.segment_size;
  int64_t pi_sqrtz = pi[isqrt(z)];
  int64_t low1 = max(low, 1);
  int64_t limit = min(low + segments * segment_size, xz);
  int64_t max_b = pi[min3(isqrt(x / low1), isqrt(limit), x_star)];
  int64_t min_b = pi[min(xz / limit, x_star)];
  min_b = max(k, min_b) + 1;

  if (min_b > max_b)
    return 0;

  Sieve sieve(low, segment_size, max_b);
  auto phi = generate_phi(phiCache, low, max_b);
  thread.init_finished();

  // Segmented sieve of Eratosthenes
  for (; low < limit; low += segment_size)
  {
    // current segment [low, high[
    int64_t high = min(low + segment_size, limit);
    low1 = max(low, 1);

    // For i < min_b there are no special leaves:
    // low <= x / (primes[i] * m) < high
    sieve.pre_sieve(primes, min_b - 1, low, high);
    int64_t b = min_b;

    // For k + 1 <= b <= pi_sqrtz
    // Find all special leaves in the current segment that are
    // composed of a prime and a square free number:
    // low <= x / (primes[b] * m) < high
    for (int64_t end = min(pi_sqrtz, max_b); b <= end; b++)
    {
      int64_t prime = primes[b];
      T xp = x / prime;
      int64_t xp_low = min(fast_div(xp, low1), z);
      int64_t xp_high = min(fast_div(xp, high), z);
      int64_t min_m = max(xp_high, z / prime);
      int64_t max_m = min(fast_div(xp, prime * prime), xp_low);

      if (prime >= max_m)
        goto next_segment;

      min_m = factor.to_index(min_m);
      max_m = factor.to_index(max_m);

      for (int64_t m = max_m; m > min_m; m--)
      {
        // mu[m] != 0 && 
        // lpf[m] > prime &&
        // mpf[m] <= y
        if (prime < factor.is_leaf(m))
        {
          int64_t xpm = fast_div64(xp, factor.to_number(m));
          int64_t stop = xpm - low;
          int64_t phi_xpm = phi[b] + sieve.count(stop);
          int64_t mu_m = factor.mu(m);
          sum -= mu_m * phi_xpm;
        }
      }

      phi[b] += sieve.get_total_count();
      sieve.cross_off_count(prime, b);
    }

    // For pi_sqrtz < b <= pi_x_star
    // Find all special leaves in the current segment
    // that are composed of 2 primes:
    // low <= x / (primes[b] * primes[l]) < high
    for (; b <= max_b; b++)
    {
      int64_t prime = primes[b];
      T xp = x / prime;
      int64_t xp_low = min(fast_div(xp, low1), y);
      int64_t xp_high = min(fast_div(xp, high), y);
      int64_t min_m = max(xp_high, prime);
      int64_t max_m = min(fast_div(xp, prime * prime), xp_low);
      int64_t l = pi[max_m];

      if (prime >= primes[l])
        goto next_segment;

      for (; primes[l] > min_m; l--)
      {
        int64_t xpq = fast_div64(xp, primes[l]);
        int64_t stop = xpq - low;
        int64_t phi_xpq = phi[b] + sieve.count(stop);
        sum += phi_xpq;
      }

      phi[b] += sieve.get_total_count();
      sieve.cross_off_count(prime, b);
    }

    next_segment:;
  }

  return sum;
}

/// D MPI worker process.
/// Asks MPI main process for new work and reports
/// partial results to MPI main process.
///
template <typename T, typename DFactorTable, typename Primes>
void D_mpi_worker(T x,
                  int64_t y,
                  int64_t z,
                  int64_t k,
                  const Primes& primes,
                  const DFactorTable& factor,
                  int threads)
{
  PiTable pi(y, threads);
  int64_t xz = x / z;
  int64_t x_star = get_x_star_gourdon(x, y);
  threads = ideal_num_threads(threads, xz);

  MpiMsg msg;
  int main_proc_id = mpi_main_proc_id();
  int proc_id = mpi_proc_id();

  #pragma omp parallel for num_threads(threads)
  for (int i = 0; i < threads; i++)
  {
    ThreadSettings thread;
    PhiCache<Primes> phiCache(primes, pi);

    while (true)
    {
      #pragma omp critical (mpi_sync)
      {
        // send result to main process
        msg.set(proc_id, i, thread.low, thread.segments, thread.segment_size, thread.sum, thread.init_secs, thread.secs);
        msg.send(main_proc_id);

        // receive new work todo
        msg.recv(proc_id);
        thread.low = msg.low();
        thread.segments = msg.segments();
        thread.segment_size = msg.segment_size();
      }

      if (thread.low >= xz)
        break;

      // Unsigned integer division is usually slightly
      // faster than signed integer division
      using UT = typename make_unsigned<T>::type;

      thread.start_time();
      UT sum = D_thread((UT) x, x_star, xz, y, z, k, primes, pi, factor, phiCache, thread);
      thread.sum = (T) sum;
      thread.stop_time();
    }
  }

  msg.set_finished();
  msg.send(main_proc_id);
}

/// D MPI main process.
/// Assigns work to the MPI worker processes.
///
template <typename T>
T D_mpi_main(T x,
             int64_t z,
             T d_approx)
{
  T sum = 0;
  int64_t xz = x / z;
  int workers = mpi_num_procs() - 1;

  MpiMsg msg;
  MpiLoadBalancer loadBalancer(x, xz, d_approx);
  Status status(x);

  while (workers > 0)
  {
    // wait for results from worker process
    msg.recv_any();

    if (msg.finished())
      workers--;
    else
    {
      sum += (T) msg.sum();
      int64_t high = msg.low() + msg.segments() * msg.segment_size();

      // update msg with new work
      loadBalancer.get_work(&msg);

      // send new work to worker process
      msg.send(msg.proc_id());
      status.print(high, xz, sum, d_approx);
    }
  }

  return sum;
}

} // namespace

namespace primecount {

int64_t D_mpi(int64_t x,
              int64_t y,
              int64_t z,
              int64_t k,
              int64_t d_approx,
              int threads)
{
  print("");
  print("=== D_mpi(x, y) ===");
  print_gourdon_vars(x, y, z, k, threads);

  int64_t sum = 0;
  double time = get_time();

  if (is_mpi_main_proc())
    sum = D_mpi_main(x, z, d_approx);
  else
  {
    DFactorTable<uint16_t> factor(y, z, threads);
    auto primes = generate_primes<int32_t>(y);
    D_mpi_worker(x, y, z, k, primes, factor, threads);
  }

  print("D", sum, time);
  return sum;
}

#ifdef HAVE_INT128_T

int128_t D_mpi(int128_t x,
               int64_t y,
               int64_t z,
               int64_t k,
               int128_t d_approx,
               int threads)
{
  print("");
  print("=== D_mpi(x, y) ===");
  print_gourdon_vars(x, y, z, k, threads);

  int128_t sum = 0;
  double time = get_time();

  if (is_mpi_main_proc())
    sum = D_mpi_main(x, z, d_approx);
  else
  {
    // uses less memory
    if (y <= FactorTable<uint16_t>::max())
    {
      DFactorTable<uint16_t> factor(y, z, threads);
      auto primes = generate_primes<uint32_t>(y);
      D_mpi_worker(x, y, z, k, primes, factor, threads);
    }
    else
    {
      DFactorTable<uint32_t> factor(y, z, threads);
      auto primes = generate_primes<int64_t>(y);
      D_mpi_worker(x, y, z, k, primes, factor, threads);
    }
  }

  print("D", sum, time);
  return sum;
}

#endif


} // namespace
