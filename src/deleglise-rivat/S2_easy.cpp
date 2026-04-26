///
/// @file  S2_easy.cpp
/// @brief Calculate the contribution of the clustered easy leaves
///        and the sparse easy leaves in parallel using OpenMP
///        (Deleglise-Rivat algorithm).
///
///        This implementation is based on the paper:
///        Tomás Oliveira e Silva, Computing pi(x): the combinatorial
///        method, Revista do DETUA, vol. 4, no. 6, March 2006,
///        pp. 759-768.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <primecount-internal.hpp>
#include <fast_div.hpp>
#include <generate_primes.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <imath.hpp>
#include <print.hpp>
#include <RelaxedAtomic.hpp>
#include <StatusS2.hpp>
#include <S.hpp>
#include <Vector.hpp>

#include <stdint.h>

#if defined(ENABLE_LIBDIVIDE)
  #include <libdivide.h>
#endif

using namespace primecount;

namespace {

#if !defined(ENABLE_LIBDIVIDE)

/// Calculate the contribution of the clustered easy leaves
/// and the sparse easy leaves.
/// @param T  either int64_t or uint128_t.
///
template <typename T, typename Primes>
T S2_easy_OpenMP(T x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 const Primes& primes,
                 int threads,
                 bool is_print)
{
  T sum = 0;
  int64_t x13 = iroot<3>(x);

  // These load balancing settings work well on my
  // dual-socket AMD EPYC 7642 server with 192 CPU cores.
  int64_t thread_threshold = 1000;
  int max_threads = (int) std::pow(z, 1 / 4.0);
  threads = std::min(threads, max_threads);
  threads = ideal_num_threads(x13, threads, thread_threshold);

  StatusS2 status(x, y, is_print);
  PiTable pi(y, threads);
  int64_t pi_sqrty = pi[isqrt(y)];
  int64_t pi_x13 = pi[x13];
  RelaxedAtomic<int64_t> min_b(max(c, pi_sqrty) + 1);

  // for (b = pi[sqrty] + 1; b <= pi_x13; b++)
  #pragma omp parallel num_threads(threads) reduction(+: sum)
  for (int64_t b = min_b++; b <= pi_x13; b = min_b++)
  {
    int64_t prime = primes[b];
    T xp = x / prime;
    int64_t min_trivial = min(xp / prime, y);
    int64_t min_clustered = (int64_t) isqrt(xp);
    int64_t min_sparse = z / prime;

    min_clustered = in_between(prime, min_clustered, y);
    min_sparse = in_between(prime, min_sparse, y);

    int64_t i = pi[min_trivial];
    int64_t pi_min_clustered = pi[min_clustered];
    int64_t pi_min_sparse = pi[min_sparse];

    if (i > pi_min_clustered)
    {
      int64_t ihi = i;
      int64_t ilo = pi_min_clustered + 1;

      // Find all clustered easy leaves where
      // successive leaves are identical.
      // pq = primes[b] * primes[i]
      // Which satisfy: pq > z && x / pq <= y
      // where phi(x / pq, b - 1) = pi(x / pq) - b + 2
      //
      // The clustered easy leaves algorithm has poor instruction
      // level parallelism because of long instruction dependency
      // chains. To mitigate this issue we process clustered
      // easy leaves bidirectionally (high-end and low-end streams)
      // which increases the number of independent instructions
      // and improves performance on modern out-of-order CPUs.
      //
      while (ilo <= ihi)
      {
        // High-end stream (decreasing i)
        int64_t xpq_hi = fast_div64(xp, primes[ihi]);
        int64_t pi_xpq_hi = pi[xpq_hi];
        int64_t phi_xpq_hi = pi_xpq_hi - b + 2;
        int64_t xpq2_hi = fast_div64(xp, primes[pi_xpq_hi + 1]);
        int64_t ihi_min = pi[xpq2_hi];
        ASSERT(ihi_min + 1 >= ilo);
        sum += phi_xpq_hi * (ihi - ihi_min);
        ihi = ihi_min;

        if (ilo > ihi)
          break;

        // Low-end stream (increasing i)
        int64_t xpq_lo = fast_div64(xp, primes[ilo]);
        int64_t pi_xpq_lo = pi[xpq_lo];
        int64_t phi_xpq_lo = pi_xpq_lo - b + 2;
        int64_t xpq2_lo = fast_div64(xp, primes[pi_xpq_lo]);
        int64_t ilo_max = pi[xpq2_lo] + 1;
        ASSERT(ilo_max - 1 <= ihi);
        sum += phi_xpq_lo * (ilo_max - ilo);
        ilo = ilo_max;
      }

      i = pi_min_clustered;
    }

    // Find all sparse easy leaves where
    // successive leaves are different.
    // pq = primes[b] * primes[i]
    // Which satisfy: pq > z && x / pq <= y
    // where phi(x / pq, b - 1) = pi(x / pq) - b + 2
    for (; i > pi_min_sparse; i--)
    {
      int64_t xpq = fast_div64(xp, primes[i]);
      sum += pi[xpq] - b + 2;
    }

    #if defined(_OPENMP) && _OPENMP >= 202011
      #pragma omp masked
    #else
      #pragma omp master
    #endif
    if (is_print)
      status.print(b, pi_x13);
  }

  return sum;
}

#elif defined(ENABLE_LIBDIVIDE)

/// This is an optimized version of S2_easy(x, y) using libdivide.
/// libdivide allows to replace expensive integer divsion
/// instructions by a sequence of shift, add and multiply
/// instructions that will calculate the integer division much
/// faster, especially on older CPUs.

/// xp < 2^64
template <typename T,
          typename LibdividePrimes>
T S2_easy_64(T xp128,
             uint64_t y,
             uint64_t z,
             uint64_t b,
             uint64_t prime,
             const LibdividePrimes& primes,
             const PiTable& pi)
{
  uint64_t xp = (uint64_t) xp128;
  uint64_t min_trivial = min(xp / prime, y);
  uint64_t min_clustered = isqrt(xp);
  uint64_t min_sparse = z / prime;
  min_clustered = in_between(prime, min_clustered, y);
  min_sparse = in_between(prime, min_sparse, y);
  uint64_t i = pi[min_trivial];
  uint64_t pi_min_clustered = pi[min_clustered];
  uint64_t pi_min_sparse = pi[min_sparse];

  T sum = 0;

  if (i > pi_min_clustered)
  {
    uint64_t ihi = i;
    uint64_t ilo = pi_min_clustered + 1;

    // Find all clustered easy leaves where
    // successive leaves are identical.
    // pq = primes[b] * primes[i]
    // Which satisfy: pq > z && x / pq <= y
    // where phi(x / pq, b - 1) = pi(x / pq) - b + 2
    //
    // The clustered easy leaves algorithm has poor instruction
    // level parallelism because of long instruction dependency
    // chains. To mitigate this issue we process clustered
    // easy leaves bidirectionally (high-end and low-end streams)
    // which increases the number of independent instructions
    // and improves performance on modern out-of-order CPUs.
    //
    while (ilo <= ihi)
    {
      // High-end stream (decreasing i)
      uint64_t xpq_hi = xp / primes[ihi];
      uint64_t pi_xpq_hi = pi[xpq_hi];
      uint64_t phi_xpq_hi = pi_xpq_hi - b + 2;
      uint64_t xpq2_hi = xp / primes[pi_xpq_hi + 1];
      uint64_t ihi_min = pi[xpq2_hi];
      ASSERT(ihi_min + 1 >= ilo);
      sum += phi_xpq_hi * (ihi - ihi_min);
      ihi = ihi_min;

      if (ilo > ihi)
        break;

      // Low-end stream (increasing i)
      uint64_t xpq_lo = xp / primes[ilo];
      uint64_t pi_xpq_lo = pi[xpq_lo];
      uint64_t phi_xpq_lo = pi_xpq_lo - b + 2;
      uint64_t xpq2_lo = xp / primes[pi_xpq_lo];
      uint64_t ilo_max = pi[xpq2_lo] + 1;
      ASSERT(ilo_max - 1 <= ihi);
      sum += phi_xpq_lo * (ilo_max - ilo);
      ilo = ilo_max;
    }

    i = pi_min_clustered;
  }

  // Find all sparse easy leaves where
  // successive leaves are different.
  // pq = primes[b] * primes[i]
  // Which satisfy: pq > z && x / pq <= y
  // where phi(x / pq, b - 1) = pi(x / pq) - b + 2
  for (; i > pi_min_sparse; i--)
  {
    uint64_t xpq = xp / primes[i];
    sum += pi[xpq] - b + 2;
  }

  return sum;
}

/// xp >= 2^64
template <typename T,
          typename Primes>
T S2_easy_128(T xp,
              uint64_t y,
              uint64_t z,
              uint64_t b,
              uint64_t prime,
              const Primes& primes,
              const PiTable& pi)
{
  uint64_t min_trivial = min(xp / prime, y);
  uint64_t min_clustered = (uint64_t) isqrt(xp);
  uint64_t min_sparse = z / prime;
  min_clustered = in_between(prime, min_clustered, y);
  min_sparse = in_between(prime, min_sparse, y);
  uint64_t i = pi[min_trivial];
  uint64_t pi_min_clustered = pi[min_clustered];
  uint64_t pi_min_sparse = pi[min_sparse];

  T sum = 0;

  if (i > pi_min_clustered)
  {
    uint64_t ihi = i;
    uint64_t ilo = pi_min_clustered + 1;

    // Find all clustered easy leaves where
    // successive leaves are identical.
    // pq = primes[b] * primes[i]
    // Which satisfy: pq > z && x / pq <= y
    // where phi(x / pq, b - 1) = pi(x / pq) - b + 2
    //
    // The clustered easy leaves algorithm has poor instruction
    // level parallelism because of long instruction dependency
    // chains. To mitigate this issue we process clustered
    // easy leaves bidirectionally (high-end and low-end streams)
    // which increases the number of independent instructions
    // and improves performance on modern out-of-order CPUs.
    //
    while (ilo <= ihi)
    {
      // High-end stream (decreasing i)
      uint64_t xpq_hi = fast_div64(xp, primes[ihi]);
      uint64_t pi_xpq_hi = pi[xpq_hi];
      uint64_t phi_xpq_hi = pi_xpq_hi - b + 2;
      uint64_t xpq2_hi = fast_div64(xp, primes[pi_xpq_hi + 1]);
      uint64_t ihi_min = pi[xpq2_hi];
      ASSERT(ihi_min + 1 >= ilo);
      sum += phi_xpq_hi * (ihi - ihi_min);
      ihi = ihi_min;

      if (ilo > ihi)
        break;

      // Low-end stream (increasing i)
      uint64_t xpq_lo = fast_div64(xp, primes[ilo]);
      uint64_t pi_xpq_lo = pi[xpq_lo];
      uint64_t phi_xpq_lo = pi_xpq_lo - b + 2;
      uint64_t xpq2_lo = fast_div64(xp, primes[pi_xpq_lo]);
      uint64_t ilo_max = pi[xpq2_lo] + 1;
      ASSERT(ilo_max - 1 <= ihi);
      sum += phi_xpq_lo * (ilo_max - ilo);
      ilo = ilo_max;
    }

    i = pi_min_clustered;
  }

  // Find all sparse easy leaves where
  // successive leaves are different.
  // pq = primes[b] * primes[i]
  // Which satisfy: pq > z && x / pq <= y
  // where phi(x / pq, b - 1) = pi(x / pq) - b + 2
  for (; i > pi_min_sparse; i--)
  {
    uint64_t xpq = fast_div64(xp, primes[i]);
    sum += pi[xpq] - b + 2;
  }

  return sum;
}

/// Calculate the contribution of the clustered easy
/// leaves and the sparse easy leaves.
///
template <typename T,
          typename Primes>
T S2_easy_OpenMP(T x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 const Primes& primes,
                 int threads,
                 bool is_print)
{
  // Initialize libdivide vector from primes vector
  Vector<libdivide::branchfree_divider<uint64_t>> lprimes;
  lprimes.resize(primes.size());
  for (std::size_t i = 1; i < lprimes.size(); i++)
    lprimes[i] = primes[i];

  T sum = 0;
  int64_t x13 = iroot<3>(x);

  // These load balancing settings work well on my
  // dual-socket AMD EPYC 7642 server with 192 CPU cores.
  int64_t thread_threshold = 1000;
  int max_threads = (int) std::pow(z, 1 / 4.0);
  threads = std::min(threads, max_threads);
  threads = ideal_num_threads(x13, threads, thread_threshold);

  StatusS2 status(x, y, is_print);
  PiTable pi(y, threads);
  int64_t pi_sqrty = pi[isqrt(y)];
  int64_t pi_x13 = pi[x13];
  RelaxedAtomic<int64_t> min_b(max(c, pi_sqrty) + 1);

  // for (b = pi[sqrty] + 1; b <= pi_x13; b++)
  #pragma omp parallel num_threads(threads) reduction(+: sum)
  for (int64_t b = min_b++; b <= pi_x13; b = min_b++)
  {
    int64_t prime = primes[b];
    T xp = x / prime;

    if (xp <= pstd::numeric_limits<uint64_t>::max())
      sum += S2_easy_64(xp, y, z, b, prime, lprimes, pi);
    else
      sum += S2_easy_128(xp, y, z, b, prime, primes, pi);

    #if defined(_OPENMP) && _OPENMP >= 202011
      #pragma omp masked
    #else
      #pragma omp master
    #endif
    if (is_print)
      status.print(b, pi_x13);
  }

  return sum;
}

#endif

void print_algo_name()
{
  #if defined(ENABLE_LIBDIVIDE)
    print("Algorithm: libdivide");
  #else
    print("Algorithm: CPU div");
  #endif
}

} // namespace

namespace primecount {

int64_t S2_easy(int64_t x,
                int64_t y,
                int64_t z,
                int64_t c,
                int threads,
                bool is_print)
{
  double time;

  if (is_print)
  {
    print("");
    print("=== S2_easy(x, y) ===");
    print_algo_name();
    print_vars(x, y, c, threads);
    time = get_time();
  }

  auto primes = generate_primes<uint32_t>(y);
  int64_t sum = S2_easy_OpenMP((uint64_t) x, y, z, c, primes, threads, is_print);

  if (is_print)
    print("S2_easy", sum, time);

  return sum;
}

#ifdef HAVE_INT128_T

int128_t S2_easy(int128_t x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 int threads,
                 bool is_print)
{
  double time;

  if (is_print)
  {
    print("");
    print("=== S2_easy(x, y) ===");
    print_algo_name();
    print_vars(x, y, c, threads);
    time = get_time();
  }

  int128_t sum;

  // uses less memory
  if (y <= pstd::numeric_limits<uint32_t>::max())
  {
    auto primes = generate_primes<uint32_t>(y);
    sum = S2_easy_OpenMP((uint128_t) x, y, z, c, primes, threads, is_print);
  }
  else
  {
    auto primes = generate_primes<int64_t>(y);
    sum = S2_easy_OpenMP((uint128_t) x, y, z, c, primes, threads, is_print);
  }

  if (is_print)
    print("S2_easy", sum, time);

  return sum;
}

#endif

} // namespace
