///
/// @file  D_arm_neon.hpp
/// @brief ARM NEON implementation of the D formula (hard special
///        leaves) in Xavier Gourdon's prime counting algorithm.
///        This algorithm is identical to D_thread_default() in D_default.hpp
///        except for NEON vectorization of index filtering.
///
///        In-depth description of this algorithm:
///        https://github.com/kimwalisch/primecount/blob/master/doc/Hard-Special-Leaves-SIMD-Filtering.pdf
///        https://github.com/kimwalisch/primecount/blob/master/doc/Hard-Special-Leaves.pdf
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef D_ARM_NEON_HPP
#define D_ARM_NEON_HPP

#include <arm_neon.h>

namespace {

using namespace primecount;

template <typename Index, int Bits>
struct NeonCompactOffsets
{
  // Selected offsets are packed first; unused slots remain zero.
  alignas(16) Index values[1 << Bits][Bits];

  constexpr NeonCompactOffsets() : values{}
  {
    for (int mask = 0; mask < (1 << Bits); mask++)
    {
      int count = 0;

      for (int bit = 0; bit < Bits; bit++)
        if (mask & (1 << bit))
          values[mask][count++] = bit;
    }
  }
};

constexpr NeonCompactOffsets<uint32_t, 8> neon_compact_offsets;
constexpr NeonCompactOffsets<uint64_t, 4> neon_compact_offsets64;

ALWAYS_INLINE uint32x4_t load_factor_u32_arm_neon(const uint16_t* factor_table)
{
  return vmovl_u16(vld1_u16(factor_table));
}

ALWAYS_INLINE uint32x4_t load_factor_u32_arm_neon(const uint32_t* factor_table)
{
  return vld1q_u32(factor_table);
}

ALWAYS_INLINE void load_factor8_u32_arm_neon(const uint16_t* factor_table,
                                            uint32x4_t& low,
                                            uint32x4_t& high)
{
  uint16x8_t factors = vld1q_u16(factor_table);
  low = vmovl_u16(vget_low_u16(factors));
  high = vmovl_high_u16(factors);
}

ALWAYS_INLINE void load_factor8_u32_arm_neon(const uint32_t* factor_table,
                                            uint32x4_t& low,
                                            uint32x4_t& high)
{
  low = vld1q_u32(factor_table);
  high = vld1q_u32(factor_table + 4);
}

/// Compute the contribution of the hard special leaves using
/// a segmented sieve. Each thread processes the interval
/// [low, low + segment_size * segments[.
///
template <typename T, typename Primes, typename FactorTable>
T D_thread_default(T x,
                   int64_t x_star,
                   int64_t xz,
                   int64_t y,
                   int64_t z,
                   int64_t k,
                   const Primes& primes,
                   const PiTable& pi,
                   const FactorTable& factor,
                   ThreadData& thread)
{
  T sum = 0;

  int64_t low = thread.low;
  int64_t low1 = max(low, 1);
  int64_t segments = thread.segments;
  int64_t segment_size = thread.segment_size;
  int64_t pi_sqrtz = pi[isqrt(z)];
  int64_t limit = min(low + segment_size * segments, xz);
  int64_t max_b = pi[min3(isqrt(x / low1), isqrt(limit), x_star)];
  int64_t min_b = pi[min(xz / limit, x_star)];
  min_b = max(k, min_b) + 1;

  if (min_b > max_b)
    return 0;

  Vector<int64_t> phi = phi_vector(low, max_b, primes, pi);
  Sieve sieve(low, segment_size, max_b);
  thread.init_time = get_time();

  INDETERMINATE Array<uint32_t, 128> m_indexes32;
  INDETERMINATE Array< int64_t, 128> m_indexes64;
  INDETERMINATE Array< int64_t, 128> xpm_low;
  const auto* factor_table = factor.data();

  // Segmented sieve of Eratosthenes
  for (; low < limit; low += segment_size)
  {
    // current segment [low, high[
    int64_t high = min(low + segment_size, limit);
    low1 = max(low, 1);

    // For b < min_b there are no special leaves:
    // low <= x / (primes[b] * m) < high
    sieve.pre_sieve(primes, min_b - 1, low, high);
    sieve.init_counter(low, high);
    int64_t b = min_b;

    // For k + 1 <= b <= pi_sqrtz
    // Find all special leaves in the current segment that are
    // composed of a prime and a square free number:
    // low <= x / (primes[b] * m) < high
    for (int64_t last = min(pi_sqrtz, max_b); b <= last; b++)
    {
      int64_t prime = primes[b];
      T xp = x / prime;
      int64_t xp_low = min(fast_div(xp, low1), z);
      int64_t xp_high = min(fast_div(xp, high), z);
      int64_t min_m = max(xp_high, z / prime);
      int64_t max_m = min(fast_div(xp, prime * prime), xp_low);

      if (prime >= max_m)
        goto next_segment;

      min_m = FactorTable::to_index(min_m);
      max_m = FactorTable::to_index(max_m);
      int64_t encoded_prime = FactorTable::encode(b);
      int64_t m = max_m;
      std::size_t m_count = 0;

      // 32-bit code path
      if (max_m <= UINT32_MAX ||
          sizeof(T) <= sizeof(uint64_t))
      {
        constexpr std::size_t max_m_count = m_indexes32.size() - 8;
        ASSERT(encoded_prime <= UINT32_MAX);
        uint32x4_t encoded_prime_vec = vdupq_n_u32(uint32_t(encoded_prime));
        // Each match contributes a mask bit and 256 to the count.
        const uint32x4_t low_weights = { 384, 320, 288, 272 };
        const uint32x4_t high_weights = { 264, 260, 258, 257 };
        uint32x4_t m_base = vdupq_n_u32(uint32_t(m));
        uint32x4_t lane_step = vdupq_n_u32(8);

        // Filter out square free m values branchlessly
        // that satisfy: factor_table[m] > encoded_prime
        for (; m >= min_m + 8; m -= 8)
        {
          uint32x4_t factor_low;
          uint32x4_t factor_high;
          load_factor8_u32_arm_neon(&factor_table[m - 7], factor_low, factor_high);
          uint32x4_t cmp_low = vcgtq_u32(factor_low, encoded_prime_vec);
          uint32x4_t cmp_high = vcgtq_u32(factor_high, encoded_prime_vec);
          uint32x4_t weighted_low = vandq_u32(cmp_low, low_weights);
          uint32x4_t weighted_high = vandq_u32(cmp_high, high_weights);
          uint32_t mask_count = vaddvq_u32(vpaddq_u32(weighted_low, weighted_high));
          const uint32_t* offsets = neon_compact_offsets.values[mask_count & 255];
          vst1q_u32(&m_indexes32[m_count], vsubq_u32(m_base, vld1q_u32(offsets)));
          vst1q_u32(&m_indexes32[m_count + 4], vsubq_u32(m_base, vld1q_u32(offsets + 4)));
          m_count += mask_count >> 8;
          m_base = vsubq_u32(m_base, lane_step);

          if (m_count > max_m_count)
          {
            // Batch calculate (xp/m - low) to improve CPU pipelining
            for (std::size_t i = 0; i < m_count; i++)
            {
              int64_t m = factor.to_number(m_indexes32[i]);
              xpm_low[i] = fast_div64(xp, m) - low;
            }

            // Process the next few special leaves that are
            // composed of a prime and a square free number:
            // low <= x / (primes[b] * m) < high
            for (std::size_t i = 0; i < m_count; i++)
            {
              // sieve.count(xp/m - low)
              int64_t count = sieve.count(xpm_low[i]);
              int64_t phi_xpm = phi[b] + count;
              sum -= factor.mu(m_indexes32[i]) * phi_xpm;
            }

            m_count = 0;
          }
        }

        // Filter out the last few square free m
        for (; m > min_m; m--)
        {
          m_indexes32[m_count] = uint32_t(m);
          m_count += (factor_table[m] > encoded_prime);
        }

        // Batch calculate (xp/m - low) to improve CPU pipelining
        for (std::size_t i = 0; i < m_count; i++)
        {
          int64_t m = factor.to_number(m_indexes32[i]);
          xpm_low[i] = fast_div64(xp, m) - low;
        }

        // Process the last few m values
        for (std::size_t i = 0; i < m_count; i++)
        {
          // sieve.count(xp/m - low)
          int64_t count = sieve.count(xpm_low[i]);
          int64_t phi_xpm = phi[b] + count;
          sum -= factor.mu(m_indexes32[i]) * phi_xpm;
        }
      }
      else // 64-bit code path
      {
        constexpr std::size_t max_m_count = m_indexes64.size() - 4;
        uint64x2_t encoded_prime_vec = vdupq_n_u64(uint64_t(encoded_prime));
        // Each match contributes a mask bit and 16 to the count.
        const uint64x2_t low_weights = { 24, 20 };
        const uint64x2_t high_weights = { 18, 17 };
        uint64x2_t m_base = vdupq_n_u64(uint64_t(m));
        uint64x2_t lane_step = vdupq_n_u64(4);

        // Filter out square free m values branchlessly
        // that satisfy: factor_table[m] > encoded_prime
        for (; m >= min_m + 4; m -= 4)
        {
          uint32x4_t factor_vec = load_factor_u32_arm_neon(&factor_table[m - 3]);
          uint64x2_t factor_low = vmovl_u32(vget_low_u32(factor_vec));
          uint64x2_t factor_high = vmovl_high_u32(factor_vec);
          uint64x2_t cmp_low = vcgtq_u64(factor_low, encoded_prime_vec);
          uint64x2_t cmp_high = vcgtq_u64(factor_high, encoded_prime_vec);
          uint64x2_t weighted_low = vandq_u64(cmp_low, low_weights);
          uint64x2_t weighted_high = vandq_u64(cmp_high, high_weights);
          uint64_t mask_count = vaddvq_u64(vaddq_u64(weighted_low, weighted_high));
          const uint64_t* offsets = neon_compact_offsets64.values[mask_count & 15];
          uint64x2_t m_low = vsubq_u64(m_base, vld1q_u64(offsets));
          uint64x2_t m_high = vsubq_u64(m_base, vld1q_u64(offsets + 2));
          vst1q_s64(&m_indexes64[m_count], vreinterpretq_s64_u64(m_low));
          vst1q_s64(&m_indexes64[m_count + 2], vreinterpretq_s64_u64(m_high));
          m_count += mask_count >> 4;
          m_base = vsubq_u64(m_base, lane_step);

          if (m_count > max_m_count)
          {
            // Batch calculate (xp/m - low) to improve CPU pipelining
            for (std::size_t i = 0; i < m_count; i++)
            {
              int64_t m = factor.to_number(m_indexes64[i]);
              xpm_low[i] = fast_div64(xp, m) - low;
            }

            // Process the next few special leaves that are
            // composed of a prime and a square free number:
            // low <= x / (primes[b] * m) < high
            for (std::size_t i = 0; i < m_count; i++)
            {
              // sieve.count(xp/m - low)
              int64_t count = sieve.count(xpm_low[i]);
              int64_t phi_xpm = phi[b] + count;
              sum -= factor.mu(m_indexes64[i]) * phi_xpm;
            }

            m_count = 0;
          }
        }

        // Filter out the last few square free m
        for (; m > min_m; m--)
        {
          m_indexes64[m_count] = m;
          m_count += (factor_table[m] > encoded_prime);
        }

        // Batch calculate (xp/m - low) to improve CPU pipelining
        for (std::size_t i = 0; i < m_count; i++)
        {
          int64_t m = factor.to_number(m_indexes64[i]);
          xpm_low[i] = fast_div64(xp, m) - low;
        }

        // Process the last few m values
        for (std::size_t i = 0; i < m_count; i++)
        {
          // sieve.count(xp/m - low)
          int64_t count = sieve.count(xpm_low[i]);
          int64_t phi_xpm = phi[b] + count;
          sum -= factor.mu(m_indexes64[i]) * phi_xpm;
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
        int64_t count = sieve.count(xpq - low);
        int64_t phi_xpq = phi[b] + count;
        sum += phi_xpq;
      }

      phi[b] += sieve.get_total_count();
      sieve.cross_off_count(prime, b);
    }

    next_segment:;
  }

  return sum;
}

} // namespace

#endif
