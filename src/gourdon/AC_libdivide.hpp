///
/// @file  AC_libdivide.hpp
/// @brief Scalar implementation of the A + C formulas in Xavier
///        Gourdon's prime counting algorithm that replaces the
///        integer divisions of numerators < 2^63 by multiplications
///        with precomputed reciprocals of the primes. This is much
///        faster than the CPU's integer division instruction,
///        especially on older CPUs. Larger numerators use the
///        shared scalar helpers in AC_default.hpp.
///
///        primecount's C formula differs from Gourdon's original
///        C formula: for the primes (x/y)^(1/3) < p <= sqrt(y),
///        the easy special leaves x/p^3 < m <= y are computed
///        together with the hard special leaves x/(p*y) < m <= x/p^3
///        using the hyperbola method, where m is a prime or a
///        product of 2 distinct primes > p. The C1 formula computes
///        the hyperbola terms at positions < p^2, which only need
///        pi[x] lookups, and the D formula computes the other terms.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef AC_LIBDIVIDE_HPP
#define AC_LIBDIVIDE_HPP

#include "AC_default.hpp"
#include "LoadBalancerAC.hpp"
#include "SegmentedPiTable.hpp"

#include <primecount-internal.hpp>
#include <fast_div.hpp>
#include <generate_primes.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <macros.hpp>
#include <min.hpp>
#include <PiTable.hpp>
#include <Vector.hpp>

#include <stdint.h>

#if defined(_MSC_VER) && \
    (defined(_M_X64) || defined(_M_ARM64)) && \
    __has_include(<intrin.h>)
  #include <intrin.h>
#endif

namespace {

using namespace primecount;

/// Returns the upper 64 bits of the 128-bit product x * y
ALWAYS_INLINE uint64_t mulhi64(uint64_t x, uint64_t y)
{
#if defined(HAVE_INT128_T)
  return uint64_t((uint128_t(x) * y) >> 64);

#elif defined(_MSC_VER) && \
      (defined(_M_X64) || defined(_M_ARM64)) && \
      __has_include(<intrin.h>)
  return __umulh(x, y);

#else
  // Portable implementation for 32-bit CPUs, computes
  // x1 * y1 * 2^64 + (x1 * y0 + x0 * y1) * 2^32 + x0 * y0.
  uint32_t x0 = uint32_t(x);
  uint32_t x1 = uint32_t(x >> 32);
  uint32_t y0 = uint32_t(y);
  uint32_t y1 = uint32_t(y >> 32);
  uint64_t x0y0 = x0 * uint64_t(y0);
  uint64_t x0y1 = x0 * uint64_t(y1);
  uint64_t x1y0 = x1 * uint64_t(y0);
  uint64_t x1y1 = x1 * uint64_t(y1);
  uint64_t mid = x1y0 + (x0y0 >> 32);
  uint64_t mid_lo = uint32_t(mid) + x0y1;
  return x1y1 + (mid >> 32) + (mid_lo >> 32);
#endif
}

/// Computes (numhi * 2^64 + numlo) / den.
/// Requires numhi < den, so that the quotient fits into 64 bits.
///
inline uint64_t div128_64(uint64_t numhi,
                          uint64_t numlo,
                          uint64_t den)
{
  ASSERT(numhi < den);

#if defined(__x86_64__) && \
   (defined(__GNUC__) || defined(__clang__))
  uint64_t quotient;
  uint64_t remainder;
  __asm__("div %[den]"
          : "=a"(quotient), "=d"(remainder)
          : [den] "r"(den), "a"(numlo), "d"(numhi));
  return quotient;

#elif defined(_MSC_VER) && \
      defined(_M_X64) && \
      _MSC_VER >= 1920 && \
      __has_include(<intrin.h>)
  uint64_t remainder;
  return _udiv128(numhi, numlo, den, &remainder);

#else
  // Knuth's Algorithm D in base 2^32, this is the same
  // algorithm as libdivide_128_div_64_to_64() in libdivide.h.
  // Normalize den so that its most significant bit is set.
  uint64_t shift = 63 - ilog2(den);
  den <<= shift;
  numhi <<= shift;
  if (shift > 0)
    numhi |= numlo >> (64 - shift);
  numlo <<= shift;

  uint64_t b = 1ull << 32;
  uint64_t den1 = den >> 32;
  uint64_t den0 = den & 0xffffffff;
  uint64_t num1 = numlo >> 32;
  uint64_t num0 = numlo & 0xffffffff;

  // Estimate the high quotient digit using den1,
  // the estimate is at most 2 too large.
  uint64_t qhat = numhi / den1;
  uint64_t rhat = numhi - qhat * den1;
  uint64_t c1 = qhat * den0;
  uint64_t c2 = rhat * b + num1;
  if (c1 > c2)
    qhat -= (c1 - c2 > den) ? 2 : 1;
  uint64_t q1 = uint32_t(qhat);

  // Estimate the low quotient digit using den1,
  // the estimate is at most 2 too large.
  uint64_t rem = numhi * b + num1 - q1 * den;
  qhat = rem / den1;
  rhat = rem - qhat * den1;
  c1 = qhat * den0;
  c2 = rhat * b + num0;
  if (c1 > c2)
    qhat -= (c1 - c2 > den) ? 2 : 1;
  uint64_t q0 = uint32_t(qhat);

  return (q1 << 32) | q0;
#endif
}

/// Divides 64-bit numerators n < 2^63 by a prime using a
/// multiplication by the prime's precomputed reciprocal,
/// based on Ben Joffe's fast floor division algorithm:
/// https://www.benjoffe.com/fast-floordiv
///
/// For d > 1, s = floor(log2(d - 1)) and m = ceil(2^(64+s) / d)
/// the error e = m * d - 2^(64+s) < d <= 2^(s+1). Hence for
/// n < 2^63, n * e < 2^(64+s) and n / d = (n * m) >> (64 + s).
///
/// Packing FastDividerU63 (9 bytes instead of 16) improves the
/// cache efficiency, this is faster even though some loads are
/// unaligned.
///
#pragma pack(push, 1)
struct FastDividerU63
{
  FastDividerU63() = default;

  FastDividerU63(uint64_t d)
  {
    ASSERT(d > 1);
    shift = uint8_t(ilog2(d - 1));
    // multiplier = ceil(2^(64+shift) / d)
    multiplier = div128_64(1ull << shift, d - 1, d);
  }

  /// Largest supported numerator
  static constexpr uint64_t max()
  {
    return pstd::numeric_limits<int64_t>::max();
  }

  uint64_t multiplier;
  uint8_t shift;
};
#pragma pack(pop)

/// Requires n <= FastDividerU63::max()
ALWAYS_INLINE uint64_t operator/(uint64_t n, const FastDividerU63& d)
{
  ASSERT(n <= FastDividerU63::max());
  return mulhi64(n, d.multiplier) >> d.shift;
}

template <typename T,
          int MULTIPLIER,
          typename LibdividePrimes>
ALWAYS_INLINE T sum_pi_libdivide(uint64_t xp,
                                 uint64_t i,
                                 uint64_t last,
                                 uint64_t b,
                                 const LibdividePrimes& primes,
                                 const SegmentedPiTable& segmentedPi)
{
  if (i > last)
    return 0;

  T sum = 0;
  uint64_t size = last - i + 1;

  // Unroll loop to increase instruction level parallelism
  for (; i + 3 <= last; i += 4)
  {
    uint64_t xpq0 = xp / primes[i];
    uint64_t xpq1 = xp / primes[i+1];
    uint64_t xpq2 = xp / primes[i+2];
    uint64_t xpq3 = xp / primes[i+3];

    sum += segmentedPi[xpq0] +
           segmentedPi[xpq1] +
           segmentedPi[xpq2] +
           segmentedPi[xpq3];
  }

  NO_UNROLL_LOOP
  for (; i <= last; i++)
  {
    uint64_t xpq = xp / primes[i];
    sum += segmentedPi[xpq];
  }

  return sum * MULTIPLIER + size * 2 - size * T(b);
}

/// Compute the A formula using FastDividerU63.
/// 64-bit numerator: xp < 2^63
/// pi[x_star] < b <= pi[x^(1/3)]
/// x / (primes[b] * primes[i]) < x^(1/2)
///
template <typename T,
          typename LibdividePrimes>
T A_libdivide(T xlow,
              T xhigh,
              uint64_t xp,
              uint64_t y,
              uint64_t prime,
              const LibdividePrimes& primes,
              const PiTable& pi,
              const SegmentedPiTable& segmentedPi)
{
  T sum = 0;

  uint64_t sqrt_xp = isqrt(xp);
  uint64_t min_2nd_prime = min(xhigh / prime, sqrt_xp);
  uint64_t max_2nd_prime = min(xlow / prime, sqrt_xp);
  uint64_t i = pi[max(prime, min_2nd_prime)] + 1;
  uint64_t max_i1 = pi[min(xp / y, max_2nd_prime)];
  uint64_t max_i2 = pi[max_2nd_prime];

  // for (; i <= max_i1; i++)
  //   sum += segmentedPi[xp / primes[i]];
  sum += sum_pi_libdivide<T, 1>(xp, i, max_i1, 2, primes, segmentedPi);
  i = max(i, max_i1 + 1);

  // for (; i <= max_i2; i++)
  //   sum += segmentedPi[xp / primes[i]] * 2;
  sum += sum_pi_libdivide<T, 2>(xp, i, max_i2, 2, primes, segmentedPi);

  return sum;
}

/// Compute the 1st part of the C formula using FastDividerU63.
/// 64-bit numerator: xp < 2^63
/// pi[(x/y)^(1/3)] < b <= pi[sqrt(y)]
/// low <= position < high
///
/// This computes the hyperbola terms at positions < p^2,
/// see C1_hyperbola() in AC_default.hpp.
///
template <typename T,
          typename LibdividePrimes>
T C1_hyperbola_libdivide(T xlow,
                         T xhigh,
                         uint64_t xp,
                         uint64_t y,
                         uint64_t b,
                         int64_t mertens_y,
                         uint64_t prime,
                         const LibdividePrimes& primes,
                         const PiTable& pi,
                         const SegmentedPiTable& segmentedPi)
{
  T sum = 0;
  uint64_t low = segmentedPi.low();
  uint64_t high = segmentedPi.high();
  uint64_t square = prime * prime;
  uint64_t xpy = xp / y;

  if (xpy >= low && xpy < high)
  {
    uint64_t phi_xpy = segmentedPi[xpy] - b + 2;
    sum -= T(mertens_y) * phi_xpy;
  }

  // If H >= p^2, only the L term has a position < p^2
  uint64_t x_div_prime3 = xp / square;

  if (x_div_prime3 >= square)
    return sum;

  uint64_t sqrt_xp = isqrt(xp);

  if (sqrt_xp >= low && sqrt_xp < high)
  {
    // M(H) * phi(H, b - 1) = (1 - u) * (u + 2)
    // = 2 - u * (u + 1) with u = pi[H] - b
    uint64_t u = segmentedPi[sqrt_xp] - b;
    sum += 2 - T(u) * (u + 1);
  }

  // Primes X/p^2 < q <= H with low <= X/q < high
  uint64_t min_q = max(fast_div64(xhigh, prime), x_div_prime3);
  uint64_t max_q = min(xlow / prime, sqrt_xp);

  if (min_q >= max_q)
    return sum;

  uint64_t min_i = pi[min_q] + 1;
  uint64_t max_i = pi[max_q];

  // for (i = min_i; i <= max_i; i++)
  //   sum += segmentedPi[xp / primes[i]] * 2 - 2 * b + 1;
  sum += sum_pi_libdivide<T, 2>(xp, min_i, max_i, 2 * b + 1, primes, segmentedPi);

  return sum;
}

/// Compute the 2nd part of the C formula using FastDividerU63.
/// Computes the clustered and sparse easy leaves of the C
/// formula for which the second factor is necessarily prime.
/// 64-bit numerator: xp < 2^63
/// pi[sqrt(y)] < b <= pi[x_star]
/// x / (primes[b] * primes[i]) < x^(1/2)
///
template <typename T,
          typename LibdividePrimes>
T C2_libdivide(T xlow,
               T xhigh,
               uint64_t xp,
               uint64_t y,
               uint64_t b,
               uint64_t pi_y,
               uint64_t max_clustered_global,
               uint64_t prime,
               const LibdividePrimes& primes,
               const PiTable& pi,
               const SegmentedPiTable& segmentedPi)
{
  uint64_t max_m = min3(xlow / prime, xp / prime, y);
  uint64_t x_div_prime3 = xp / (prime * prime);
  uint64_t xhigh_div_prime = fast_div64(xhigh, prime);
  uint64_t min_m = max3(xhigh_div_prime, x_div_prime3, prime);

  if (min_m >= max_m)
    return 0;

  uint64_t pi_min_m = pi[min_m];
  uint64_t sqrt_xp = isqrt(xp);
  uint64_t min_clustered = in_between(min_m, sqrt_xp, max_m);
  uint64_t pi_min_clustered = pi[min_clustered];
  uint64_t min_clustered_global = max3(x_div_prime3, sqrt_xp, prime);
  min_clustered_global = min(min_clustered_global, y);
  uint64_t pi_conj_lo = pi_min_m;
  uint64_t pi_conj_hi = pi_min_m;
  uint64_t i = pi_min_m + 1;

  T sum = 0;

  // Compute the boundary correction once
  if (max_m >= max_clustered_global &&
      pi_y > pi_min_clustered)
  {
    uint64_t q_lo = fast_div64(xp, max_clustered_global);
    uint64_t q_hi = fast_div64(xp, min_clustered_global + 1);
    uint64_t pi_min_clustered_global = pi[min_clustered_global];
    sum += T(pi[q_lo]) * pi_y - T(pi[q_hi]) * pi_min_clustered_global;
    sum -= T(b - 2) * (pi_y - pi_min_clustered_global);
  }

  // Reflected range: ]pi_conj_lo, pi_conj_hi]
  if (min_clustered_global < max_clustered_global &&
      segmentedPi.low() < max_clustered_global &&
      segmentedPi.high() > min_clustered_global + 1)
  {
    uint64_t q_lo = fast_div64(xp, max_clustered_global);
    uint64_t q_hi = fast_div64(xp, min_clustered_global + 1);
    pi_conj_lo = max(pi[q_lo], pi_min_m);
    pi_conj_hi = min(pi[q_hi], pi_min_clustered);
    pi_conj_hi = max(pi_conj_hi, pi_conj_lo);
  }

  // Sparse leaves below the reflected range.
  // for (; i <= pi_conj_lo; i++)
  //   sum += segmentedPi[xp / primes[i]] - b + 2;
  sum += sum_pi_libdivide<T, 1>(xp, i, pi_conj_lo, b, primes, segmentedPi);
  i = pi_conj_lo + 1;

  // Reflected leaves are counted once as a
  // sparse leaf and once as a conjugate.
  // for (; i <= pi_conj_hi; i++)
  //   sum += segmentedPi[xp / primes[i]] * 2 - b + 2;
  sum += sum_pi_libdivide<T, 2>(xp, i, pi_conj_hi, b, primes, segmentedPi);
  i = pi_conj_hi + 1;

  // Sparse leaves above the reflected range.
  // for (; i <= pi_min_clustered; i++)
  //   sum += segmentedPi[xp / primes[i]] - b + 2;
  sum += sum_pi_libdivide<T, 1>(xp, i, pi_min_clustered, b, primes, segmentedPi);

  return sum;
}

/// Compute A + C
template <typename T,
          typename Primes>
T AC_OpenMP_libdivide(T x,
                      int64_t y,
                      int64_t k,
                      int64_t x_star,
                      const PiTable& pi,
                      const Primes& primes,
                      int threads,
                      bool is_print)
{
  T sum = 0;
  int64_t x13 = iroot<3>(x);
  int64_t sqrtx = isqrt(x);
  int64_t xy = x / y;

  // These load balancing settings work well on my
  // dual-socket AMD EPYC 7642 server with 192 CPU cores.
  int64_t thread_threshold = 1000;
  int max_threads = (int) std::pow(xy, 1 / 3.7);
  threads = min(threads, max_threads);
  threads = ideal_num_threads(x13, threads, thread_threshold);
  INDETERMINATE LoadBalancerAC loadBalancer(sqrtx, y, threads, is_print);

  int64_t pi_y = pi[y];
  int64_t max_clustered_global = pi.prev_prime(y);
  int64_t sqrty = isqrt(y);
  int64_t pi_sqrty = pi[sqrty];
  int64_t pi_root3_xy = pi[iroot<3>(xy)];
  auto mertens_y = generate_mertens_y(y, k, primes, pi);

  // Initialize FastDividerU63 vector from primes vector
  Vector<FastDividerU63> lprimes;
  lprimes.resize(primes.size());

  int64_t min_thread_size = (int64_t) 1e6;
  int64_t primes_size = lprimes.size();
  int init_threads = ideal_num_threads(primes_size, threads, min_thread_size);
  int64_t thread_dist = ceil_div(primes_size, init_threads);

  #pragma omp parallel for num_threads(init_threads) schedule(static, 1)
  for (int64_t low = 1; low < primes_size; low += thread_dist)
  {
    int64_t high = min(low + thread_dist, primes_size);
    for (int64_t i = low; i < high; i++)
      lprimes[i] = FastDividerU63(primes[i]);
  }

  // In order to reduce the thread creation & destruction
  // overhead we reuse the same threads throughout the
  // entire computation. The same threads are used for:
  //
  // 1) Computation of the C1 formula.
  // 2) Computation of the C2 formula.
  // 3) Computation of the A formula.
  //
  #pragma omp parallel num_threads(threads) reduction(+: sum)
  {
    // SegmentedPiTable is accessed very frequently.
    // In order to get good performance it is important that
    // SegmentedPiTable fits into the CPU's cache.
    // Hence we use a small segment_size of x^(1/4).
    SegmentedPiTable segmentedPi;
    ThreadDataAC thread;

    // for (low = 0; low < sqrt(x); low += segment_size)
    while (loadBalancer.get_work(thread))
    {
      int64_t low = thread.low;
      int64_t segment_size = thread.segment_size;
      int64_t limit = low + thread.segments * segment_size;
      limit = min(limit, sqrtx);

      for (; low < limit; low += segment_size)
      {
        // Current segment [low, high[
        int64_t high = low + segment_size;
        high = min(high, sqrtx);
        segmentedPi.init(low, high, limit);

        // We measure the thread computation time excluding the
        // first expensive initialization of the segmentedPi
        // lookup table. If the thread computation time is close
        // to 0 then we increase the number of segments in the
        // loadBalancer which should improve performance.
        if (low == thread.low)
          thread.secs = get_time();

        int64_t pi_sqrt_low = pi[isqrt(low)];
        T xlow = x / max(low, 1);
        T xhigh = x / high;

        // The C1 hyperbola terms have positions < p^2 <= y
        if (low < y)
        {
          int64_t min_c1 = max(k, pi_root3_xy);
          int64_t min_c1_prime = min(xy / high, sqrty);
          min_c1 = max3(min_c1, pi_sqrt_low, pi[min_c1_prime]) + 1;

          // C1 formula: pi[(x/y)^(1/3)] < b <= pi[sqrt(y)]
          for (int64_t b = min_c1; b <= pi_sqrty; b++)
          {
            int64_t prime = primes[b];
            T xp = x / prime;

            if (xp <= FastDividerU63::max())
              sum += C1_hyperbola_libdivide(xlow, xhigh, uint64_t(xp), y, b, mertens_y[b], prime, lprimes, pi, segmentedPi);
            else
              sum += C1_hyperbola(xlow, xhigh, xp, y, b, mertens_y[b], primes, pi, segmentedPi);
          }
        }

        int64_t min_c2 = max3(k, pi_root3_xy, pi_sqrty);
        int64_t min_c2_prime = min(xy / high, x_star);
        min_c2 = max3(min_c2, pi_sqrt_low, pi[min_c2_prime]) + 1;
        int64_t x_div_high2 = fast_div64(xhigh, high);
        int64_t min_a = pi[in_between(x_star, x_div_high2, x13)] + 1;

        // Upper bound of A & C2 formulas:
        // x / (p * q) >= low
        // p * next_prime(p) <= x / low
        // p <= sqrt(x / low)
        int64_t sqrt_xlow = isqrt(xlow);
        int64_t max_c2 = pi[min(sqrt_xlow, x_star)];
        T max_c2_prime = xlow / max(max_clustered_global, 1);
        int64_t max_c2_clustered = pi[min3(max_c2_prime, sqrt_xlow, x_star)];
        int64_t min_c2_sparse = pi[min(x_div_high2, x_star)] + 1;
        min_c2_sparse = max3(min_c2, min_c2_sparse, max_c2_clustered + 1);
        int64_t max_a = pi[min(sqrt_xlow, x13)];

        // C2 formula: pi[sqrt(y)] < b <= pi[x_star]
        for (int64_t b = min_c2; b <= max_c2_clustered; b++)
        {
          int64_t prime = primes[b];
          T xp = x / prime;

          if (xp <= FastDividerU63::max())
            sum += C2_libdivide(xlow, xhigh, uint64_t(xp), y, b, pi_y, max_clustered_global, prime, lprimes, pi, segmentedPi);
          else
            sum += C2(xlow, xhigh, xp, y, b, pi_y, max_clustered_global, primes, pi, segmentedPi);
        }

        // C2 formula: pi[sqrt(y)] < b <= pi[x_star]
        for (int64_t b = min_c2_sparse; b <= max_c2; b++)
        {
          int64_t prime = primes[b];
          T xp = x / prime;

          if (xp <= FastDividerU63::max())
            sum += C2_libdivide(xlow, xhigh, uint64_t(xp), y, b, pi_y, max_clustered_global, prime, lprimes, pi, segmentedPi);
          else
            sum += C2(xlow, xhigh, xp, y, b, pi_y, max_clustered_global, primes, pi, segmentedPi);
        }

        // A formula: pi[x_star] < b <= pi[x13]
        for (int64_t b = min_a; b <= max_a; b++)
        {
          int64_t prime = primes[b];
          T xp = x / prime;

          if (xp <= FastDividerU63::max())
            sum += A_libdivide(xlow, xhigh, uint64_t(xp), y, prime, lprimes, pi, segmentedPi);
          else
            sum += A(xlow, xhigh, xp, y, b, primes, pi, segmentedPi);
        }
      }
    }
  }

  return sum;
}

} // namespace

#endif
