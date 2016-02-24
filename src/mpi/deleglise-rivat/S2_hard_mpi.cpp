///
/// @file  S2_hard_mpi.cpp
/// @brief Calculate the contribution of the hard special leaves which
///        require use of a sieve (Deleglise-Rivat algorithm).
///        This is a parallel implementation which uses compression
///        (PiTable & FactorTable) to reduce the memory usage by
///        about 10x.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <S2.hpp>
#include <FactorTable.hpp>
#include <primecount-internal.hpp>
#include <BitSieve.hpp>
#include <fast_div.hpp>
#include <generate.hpp>
#include <int128.hpp>
#include <min_max.hpp>
#include <pmath.hpp>
#include <S2LoadBalancer.hpp>
#include <S2Status.hpp>
#include <tos_counters.hpp>
#include <Wheel.hpp>

#include <stdint.h>
#include <vector>

#include <mpi.h>
#include <S2_hard_mpi_msg.hpp>
#include <S2_hard_mpi_LoadBalancer.hpp>
#include <mpi_reduce_sum.hpp>

#ifdef _OPENMP
  #include <omp.h>
#endif

/// Count the number of unsieved elements inside
/// [start, stop] from the sieve array.
///
#define COUNT_POPCNT(start, stop) \
    sieve.count(start, stop, low, high, count, count_low_high)

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
/// the special counters tree data structure.
///
template <typename T>
void cross_off(BitSieve& sieve,
               int64_t low,
               int64_t high,
               int64_t prime,
               WheelItem& w,
               T& counters)
{
  int64_t segment_size = sieve.size();
  int64_t m = w.next_multiple;
  int64_t wheel_index = w.wheel_index;

  for (; m < high; m += prime * Wheel::next_multiple_factor(&wheel_index))
  {
    if (sieve[m - low])
    {
      sieve.unset(m - low);
      cnt_update(counters, m - low, segment_size);
    }
  }

  w.set(m, wheel_index);
}

/// @return  true if the interval [low, high] contains
///          few hard special leaves.
///
bool is_popcnt(int64_t low,
               int64_t high,
               int64_t y,
               double alpha)
{
  return (high < y || low > y * alpha);
}

/// Compute the S2 contribution of the hard special leaves which
/// require use of a sieve. Each thread processes the interval
/// [low_thread, low_thread + segments * segment_size[
/// and the missing special leaf contributions for the interval
/// [1, low_process[ are later reconstructed and added in
/// the parent S2_hard_OpenMP_master() function.
///
template <typename T, typename FactorTable, typename Primes>
T S2_hard_OpenMP_thread(T x,
                        int64_t y,
                        int64_t z,
                        int64_t c,
                        int64_t segment_size,
                        int64_t segments_per_thread,
                        int64_t thread_num,
                        int64_t low,
                        int64_t limit,
                        double alpha,
                        FactorTable& factors,
                        PiTable& pi,
                        Primes& primes,
                        vector<int64_t>& mu_sum,
                        vector<int64_t>& phi)
{
  low += segment_size * segments_per_thread * thread_num;
  limit = min(low + segment_size * segments_per_thread, limit);
  int64_t max_b = pi[min3(isqrt(x / low), isqrt(z), y)];
  int64_t pi_sqrty = pi[isqrt(y)];
  T s2_hard = 0;

  if (c > max_b)
    return s2_hard;

  BitSieve sieve(segment_size);
  Wheel wheel(primes, max_b + 1, low);
  vector<int32_t> counters;
  phi.resize(max_b + 1, 0);
  mu_sum.resize(max_b + 1, 0);

  // Segmented sieve of Eratosthenes
  for (; low < limit; low += segment_size)
  {
    // Current segment = interval [low, high[
    int64_t high = min(low + segment_size, limit);
    int64_t b = c + 1;

    // pre-sieve the multiples of the first c primes
    sieve.pre_sieve(c, low);

    // Calculate the contribution of the hard special leaves using the
    // POPCNT algorithm. If there are relatively few special leaves
    // per segment we count the number of unsieved elements directly
    // from the sieve array using the POPCNT instruction.
    if (is_popcnt(low, high, y, alpha))
    {
      int64_t count_low_high = sieve.count((high - 1) - low);

      // For c + 1 <= b <= pi_sqrty
      // Find all special leaves: n = primes[b] * m
      // which satisfy: mu[m] != 0 && primes[b] < lpf[m] && low <= (x / n) < high
      for (int64_t end = min(pi_sqrty, max_b); b <= end; b++)
      {
        int64_t prime = primes[b];
        T x2 = x / prime;
        int64_t x2_div_low = min(fast_div(x2, low), y);
        int64_t x2_div_high = min(fast_div(x2, high), y);
        int64_t min_m = max(x2_div_high, y / prime);
        int64_t max_m = x2_div_low;
        int64_t count = 0;
        int64_t start = 0;

        if (prime >= max_m)
          goto next_segment;

        factors.to_index(&min_m);
        factors.to_index(&max_m);

        for (int64_t m = max_m; m > min_m; m--)
        {
          if (prime < factors.lpf(m))
          {
            int64_t xn = (int64_t) fast_div(x2, factors.get_number(m));
            int64_t stop = xn - low;
            count += COUNT_POPCNT(start, stop);
            start = stop + 1;
            int64_t phi_xn = phi[b] + count;
            int64_t mu_m = factors.mu(m);
            s2_hard -= mu_m * phi_xn;
            mu_sum[b] -= mu_m;
          }
        }

        phi[b] += count_low_high;
        count_low_high -= cross_off(sieve, low, high, prime, wheel[b]);
      }

      // For pi_sqrty <= b <= pi_sqrtz
      // Find all hard special leaves: n = primes[b] * primes[l]
      // which satisfy: low <= (x / n) < high
      for (; b <= max_b; b++)
      {
        int64_t prime = primes[b];
        T x2 = x / prime;
        int64_t x2_div_low = min(fast_div(x2, low), y);
        int64_t x2_div_high = min(fast_div(x2, high), y);
        int64_t l = pi[min(x2_div_low, z / prime)];
        int64_t min_hard_leaf = max3(x2_div_high, y / prime, prime);
        int64_t count = 0;
        int64_t start = 0;

        if (prime >= primes[l])
          goto next_segment;

        for (; primes[l] > min_hard_leaf; l--)
        {
          int64_t xn = (int64_t) fast_div(x2, primes[l]);
          int64_t stop = xn - low;
          count += COUNT_POPCNT(start, stop);
          start = stop + 1;
          int64_t phi_xn = phi[b] + count;
          s2_hard += phi_xn;
          mu_sum[b]++;
        }

        phi[b] += count_low_high;
        count_low_high -= cross_off(sieve, low, high, prime, wheel[b]);
      }
    }
    else
    {
      // Calculate the contribution of the hard special leaves using
      // Tomás Oliveira's O(log(N)) special tree data structure
      // for counting the number of unsieved elements. This algorithm
      // runs fastest if there are many special leaves per segment.

      // allocate memory upon first usage
      counters.resize(segment_size);

      // Initialize special tree data structure from sieve
      cnt_finit(sieve, counters, segment_size);

      // For c + 1 <= b <= pi_sqrty
      // Find all special leaves: n = primes[b] * m
      // which satisfy: mu[m] != 0 && primes[b] < lpf[m] && low <= (x / n) < high
      for (int64_t end = min(pi_sqrty, max_b); b <= end; b++)
      {
        int64_t prime = primes[b];
        T x2 = x / prime;
        int64_t x2_div_low = min(fast_div(x2, low), y);
        int64_t x2_div_high = min(fast_div(x2, high), y);
        int64_t min_m = max(x2_div_high, y / prime);
        int64_t max_m = x2_div_low;

        if (prime >= max_m)
          goto next_segment;

        factors.to_index(&min_m);
        factors.to_index(&max_m);

        for (int64_t m = max_m; m > min_m; m--)
        {
          if (prime < factors.lpf(m))
          {
            int64_t xn = (int64_t) fast_div(x2, factors.get_number(m));
            int64_t count = cnt_query(counters, xn - low);
            int64_t phi_xn = phi[b] + count;
            int64_t mu_m = factors.mu(m);
            s2_hard -= mu_m * phi_xn;
            mu_sum[b] -= mu_m;
          }
        }

        phi[b] += cnt_query(counters, (high - 1) - low);
        cross_off(sieve, low, high, prime, wheel[b], counters);
      }

      // For pi_sqrty <= b <= pi_sqrtz
      // Find all hard special leaves: n = primes[b] * primes[l]
      // which satisfy: low <= (x / n) < high
      for (; b <= max_b; b++)
      {
        int64_t prime = primes[b];
        T x2 = x / prime;
        int64_t x2_div_low = min(fast_div(x2, low), y);
        int64_t x2_div_high = min(fast_div(x2, high), y);
        int64_t l = pi[min(x2_div_low, z / prime)];
        int64_t min_hard_leaf = max3(x2_div_high, y / prime, prime);

        if (prime >= primes[l])
          goto next_segment;

        for (; primes[l] > min_hard_leaf; l--)
        {
          int64_t xn = (int64_t) fast_div(x2, primes[l]);
          int64_t count = cnt_query(counters, xn - low);
          int64_t phi_xn = phi[b] + count;
          s2_hard += phi_xn;
          mu_sum[b]++;
        }

        phi[b] += cnt_query(counters, (high - 1) - low);
        cross_off(sieve, low, high, prime, wheel[b], counters);
      }
    }

    next_segment:;
  }

  return s2_hard;
}

/// Calculate the contribution of the hard special leaves which
/// require use of a sieve (to reduce the memory usage).
/// This is a parallel implementation with advanced load balancing.
/// As most special leaves tend to be in the first segments we
/// start off with a small segment size and few segments
/// per thread, after each iteration we dynamically increase
/// the segment size and the segments per thread.
///
template <typename T, typename FactorTable, typename Primes>
T S2_hard_OpenMP_master(int64_t low,
                        int64_t high,
                        T x,
                        int64_t y,
                        int64_t z,
                        int64_t c,
                        int64_t segment_size,
                        int64_t segments_per_thread,
                        T s2_hard_approx,
                        Primes& primes,
                        PiTable& pi,
                        FactorTable& factors,
                        int proc_id,
                        int threads)
{
  double time = get_wtime();
  threads = validate_threads(threads, z);

  T s2_hard = 0;
  int64_t limit = high + 1;
  int64_t old_low = low;
  int64_t old_high = high;

  double alpha = get_alpha(x, y);
  S2Status status(x);
  S2LoadBalancer loadBalancer(x, y, z, threads);
  int64_t min_segment_size = loadBalancer.get_min_segment_size();

  int64_t max_b = pi[min3(isqrt(x / low), isqrt(z), y)];
  vector<int64_t> phi_total;// = phi_vector(low - 1, max_b, primes, pi, threads);

  while (low < limit)
  {
    // make sure we use all CPU cores
    segment_size = min(segment_size, (limit - low) / threads);
    segment_size = max(segment_size, min_segment_size);
    segment_size = prev_power_of_2(segment_size);

    int64_t segments = ceil_div(limit - low, segment_size);
    threads = in_between(1, threads, segments);
    segments_per_thread = in_between(1, segments_per_thread, ceil_div(segments, threads));

    aligned_vector<vector<int64_t> > phi(threads);
    aligned_vector<vector<int64_t> > mu_sum(threads);
    aligned_vector<double> timings(threads);

    #pragma omp parallel for num_threads(threads) reduction(+: s2_hard)
    for (int i = 0; i < threads; i++)
    {
      timings[i] = get_wtime();
      s2_hard += S2_hard_OpenMP_thread(x, y, z, c, segment_size, segments_per_thread,
          i, low, limit, alpha, factors, pi, primes, mu_sum[i], phi[i]);
      timings[i] = get_wtime() - timings[i];
    }

    // Once all threads have finished reconstruct and add the 
    // missing contribution of all special leaves. This must
    // be done in order as each thread (i) requires the sum of
    // the phi values from the previous threads.
    //
    for (int i = 0; i < threads; i++)
    {
      for (size_t j = 1; j < phi[i].size(); j++)
      {
        s2_hard += phi_total[j] * (T) mu_sum[i][j];
        phi_total[j] += phi[i][j];
      }
    }

    low += segments_per_thread * threads * segment_size;
    loadBalancer.update(low, threads, &segment_size, &segments_per_thread, timings);
  }

  S2_hard_mpi_msg result_msg(proc_id, old_low, old_high, segment_size,
      segments_per_thread, s2_hard, get_wtime() - time,
          loadBalancer.get_rsd());

  int master_proc_id = 0;
  result_msg.send(master_proc_id);

  return s2_hard;
}

/// S2_hard MPI slave process.
/// Computes a part of the hard speacial leaves on a cluster node
/// and sends the result to the master process.
///
template <typename F, typename T>
void S2_hard_mpi_slave(T x,
                       int64_t y,
                       int64_t z,
                       int64_t c,
                       T s2_hard_approx,
                       int proc_id,
                       int threads)
{
  // this will take a while to initialize
  FactorTable<F> factors(y);
  int64_t max_prime = z / isqrt(y);
  vector<int64_t> primes = generate_primes<int64_t>(max_prime);
  PiTable pi(max_prime);

  S2_hard_mpi_msg get_work;
  get_work.recv(proc_id);

  while (!get_work.finished())
  {
    proc_id = get_work.proc_id();

    S2_hard_OpenMP_master(get_work.low(), get_work.high(), x, y, z, c,
        get_work.segment_size(), get_work.segments_per_thread(), s2_hard_approx,
            primes, pi, factors, proc_id, threads);

    get_work.recv(proc_id);
  }
}

/// S2_hard MPI master process.
/// Distributes the computation of the hard speacial leaves on
/// cluster nodes.
///
template <typename T>
T S2_hard_mpi_master(T x,
                     int64_t y,
                     int64_t z,
                     int64_t c,
                     T s2_hard_approx,
                     int procs,
                     int threads)
{
  T s2_hard = 0;
  int slave_procs = procs - 1;

  int64_t low = 1;
  int64_t high = 0;
  int64_t sqrtz = isqrt(z);
  int64_t logx = max(1, ilog(x));

  // start with tiny segment size
  int64_t segment_size = max(sqrtz / logx, 1 << 9);
  int64_t segments_per_thread = 1;
  int64_t proc_interval = sqrtz;

  // start all slave processes
  for (int proc_id = 1; proc_id <= slave_procs; proc_id++)
  {
    low = high + 1;
    high = min(low + proc_interval, z);

    S2_hard_mpi_msg msg(proc_id, low, high, segment_size, segments_per_thread);
    msg.send(proc_id);
  }

  S2Status status(x);
  S2_hard_mpi_LoadBalancer loadBalancer(low, high, y, z, slave_procs);

  // main process scheduling loop
  while (true)
  {
    // wait for results from slave process
    S2_hard_mpi_msg msg;
    msg.recv_any();
    s2_hard += msg.s2_hard<T>();

    if (print_status())
      status.print(s2_hard, s2_hard_approx, msg.rsd());

    // update msg with new work to do
    loadBalancer.update(&msg, status.get_percent());

    if (loadBalancer.finished())
    {
      msg.send_finish();
      break;
    }
 
    // send new work to slave process
    msg.send(msg.proc_id());
  }

  // we are nearly finished, wait for remaining results
  for (int i = 1; i <= slave_procs - 1; i++)
  {
    S2_hard_mpi_msg msg;
    msg.recv_any();
    s2_hard += msg.s2_hard<T>();

    if (print_status())
      status.print(s2_hard, s2_hard_approx, msg.rsd());

    msg.send_finish();
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
  int proc_id;
  int procs;
  int master_proc_id = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  if (procs < 2)
    return S2_hard(x, y, z, c, s2_hard_approx, threads);

  print("");
  print("=== S2_hard_mpi(x, y) ===");
  print("Computation of the hard special leaves");
  print(x, y, c, threads);

  int64_t s2_hard = 0;
  double time = get_wtime();

  if (proc_id == master_proc_id)
    s2_hard = S2_hard_mpi_master(x, y, z, c, s2_hard_approx, procs, threads);
  else
    S2_hard_mpi_slave<uint16_t>(x, y, z, c, s2_hard_approx, proc_id, threads);

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
  int proc_id;
  int procs;
  int master_proc_id = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  if (procs < 2)
    return S2_hard(x, y, z, c, s2_hard_approx, threads);

  print("");
  print("=== S2_hard_mpi(x, y) ===");
  print("Computation of the hard special leaves");
  print(x, y, c, threads);

  int128_t s2_hard = 0;
  double time = get_wtime();

  if (proc_id == master_proc_id)
    s2_hard = S2_hard_mpi_master(x, y, z, c, s2_hard_approx, procs, threads);
  else
  {
    // uses less memory
    if (y <= FactorTable<uint16_t>::max())
      S2_hard_mpi_slave<uint16_t>((intfast128_t) x, y, z, c, (intfast128_t) s2_hard_approx, proc_id, threads);
    else
      S2_hard_mpi_slave<uint32_t>((intfast128_t) x, y, z, c, (intfast128_t) s2_hard_approx, proc_id, threads);
  }

  print("S2_hard", s2_hard, time);
  return s2_hard;
}

#endif

} // namespace primecount
