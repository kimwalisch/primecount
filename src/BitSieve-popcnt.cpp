///
/// @file  BitSieve-popcnt.cpp
/// @brief Count the number of 1 bits inside a 64-bit array.
///        The AVX2 popcount algorithm used in this file is described
///        in the paper "Faster Population Counts using AVX2
///        Instructions" by Wojciech Muła, Nathan Kurz, Daniel Lemire.
/// @see   https://arxiv.org/abs/1611.07612
/// @see   https://github.com/WojciechMula/sse-popcount
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <BitSieve.hpp>
#include <popcnt.hpp>

#include <stdint.h>
#include <cassert>
#include <vector>

namespace {
namespace POPCNT {

/// Count the number of 1 bits inside the data
/// array using the POPCNT instruction.
///
uint64_t popcnt(const uint64_t* data, uint64_t size)
{
  uint64_t sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;
  uint64_t limit = size - size % 4;
  uint64_t i = 0;

  for (; i < limit; i += 4)
  {
    sum0 += popcnt64(data[i+0]);
    sum1 += popcnt64(data[i+1]);
    sum2 += popcnt64(data[i+2]);
    sum3 += popcnt64(data[i+3]);
  }

  uint64_t total = sum0 + sum1 + sum2 + sum3;

  for (; i < size; i++)
    total += popcnt64(data[i]);

  return total;
}

} // namespace POPCNT
} // namespace

#if defined(HAVE_AVX2)

#include <immintrin.h>

namespace {
namespace AVX2 {

__attribute__ ((target ("avx2")))
__m256i popcnt(const __m256i v)
{
  __m256i m1 = _mm256_set1_epi8(0x55);
  __m256i m2 = _mm256_set1_epi8(0x33);
  __m256i m4 = _mm256_set1_epi8(0x0F);
  __m256i t1 = _mm256_sub_epi8(v, (_mm256_srli_epi16(v,  1) & m1));
  __m256i t2 = _mm256_add_epi8(t1 & m2, (_mm256_srli_epi16(t1, 2) & m2));
  __m256i t3 = _mm256_add_epi8(t2, _mm256_srli_epi16(t2, 4)) & m4;

  return _mm256_sad_epu8(t3, _mm256_setzero_si256());
}

__attribute__ ((target ("avx2")))
void CSA(__m256i& h, __m256i& l, __m256i a, __m256i b, __m256i c)
{
  __m256i u = a ^ b;
  h = (a & b) | (u & c);
  l = u ^ c;
}

/// AVX2 Harley-Seal popcount (4th iteration).
/// The algorithm is based on the paper "Faster Population Counts
/// using AVX2 Instructions" by Daniel Lemire, Nathan Kurz and
/// Wojciech Mula (23 Nov 2016).
/// @see https://arxiv.org/abs/1611.07612
///
__attribute__ ((target ("avx2")))
uint64_t popcnt(const __m256i* data, uint64_t size)
{
  __m256i total = _mm256_setzero_si256();
  __m256i ones = _mm256_setzero_si256();
  __m256i twos = _mm256_setzero_si256();
  __m256i fours = _mm256_setzero_si256();
  __m256i eights = _mm256_setzero_si256();
  __m256i sixteens = _mm256_setzero_si256();
  __m256i twosA, twosB, foursA, foursB, eightsA, eightsB;

  uint64_t limit = size - size % 16;
  uint64_t i = 0;

  for(; i < limit; i += 16)
  {
    CSA(twosA, ones, ones, data[i+0], data[i+1]);
    CSA(twosB, ones, ones, data[i+2], data[i+3]);
    CSA(foursA, twos, twos, twosA, twosB);
    CSA(twosA, ones, ones, data[i+4], data[i+5]);
    CSA(twosB, ones, ones, data[i+6], data[i+7]);
    CSA(foursB, twos, twos, twosA, twosB);
    CSA(eightsA,fours, fours, foursA, foursB);
    CSA(twosA, ones, ones, data[i+8], data[i+9]);
    CSA(twosB, ones, ones, data[i+10], data[i+11]);
    CSA(foursA, twos, twos, twosA, twosB);
    CSA(twosA, ones, ones, data[i+12], data[i+13]);
    CSA(twosB, ones, ones, data[i+14], data[i+15]);
    CSA(foursB, twos, twos, twosA, twosB);
    CSA(eightsB, fours, fours, foursA, foursB);
    CSA(sixteens, eights, eights, eightsA, eightsB);

    total = _mm256_add_epi64(total, popcnt(sixteens));
  }

  total = _mm256_slli_epi64(total, 4);
  total = _mm256_add_epi64(total, _mm256_slli_epi64(popcnt(eights), 3));
  total = _mm256_add_epi64(total, _mm256_slli_epi64(popcnt(fours), 2));
  total = _mm256_add_epi64(total, _mm256_slli_epi64(popcnt(twos), 1));
  total = _mm256_add_epi64(total, popcnt(ones));

  for(; i < size; i++)
    total = _mm256_add_epi64(total, popcnt(data[i]));

  uint64_t* total64 = (uint64_t*) &total;

  return total64[0] +
         total64[1] +
         total64[2] +
         total64[3];
}

/// Align memory to 32 bytes boundary
void align(const uint64_t*& p,
           uint64_t* size,
           uint64_t* total)
{
  for (; (uintptr_t) p % 32 != 0; p++)
  {
    assert(*size > 0);
    *total += popcnt64(*p);
    *size -= 1;
  }
}

/// AVX2 popcount algorithm.
/// @param data  A 64-bit array
/// @param size  Length of data array
///
uint64_t popcnt(const uint64_t* data, uint64_t size)
{
  uint64_t total = 0;

  // AVX2 popcount is faster than POPCNT 
  // for array sizes >= 1 kilobyte
  if (size * 8 >= 1024)
  {
    align(data, &size, &total);
    total += popcnt((const __m256i*) data, size / 4);
    data += size - size % 4;
    size = size % 4;
  }

  // process remaining words
  total += POPCNT::popcnt(data, size);

  return total;
}

} // namespace AVX2
} // namespace

#endif

namespace {

#if defined(HAVE_AVX2)

// calls CPUID at program startup
const int avx2 = __builtin_cpu_supports("avx2");

uint64_t popcnt(const uint64_t* data, uint64_t size)
{
  if (avx2)
    return AVX2::popcnt(data, size);
  else
    return POPCNT::popcnt(data, size);
}

#else

uint64_t popcnt(const uint64_t* data, uint64_t size)
{
  return POPCNT::popcnt(data, size);
}

#endif

} // namespace

namespace primecount {

/// Count the number of 1 bits inside [start, stop]
uint64_t BitSieve::count(uint64_t start,
                         uint64_t stop) const
{
  if (start > stop)
    return 0;

  assert(stop < size_);

  uint64_t start_idx = start / 64;
  uint64_t stop_idx = stop / 64;
  uint64_t m1 = 0xffffffffffffffffull << (start % 64);
  uint64_t m2 = 0xffffffffffffffffull >> (63 - stop % 64);
  uint64_t bit_count;

  if (start_idx == stop_idx)
    bit_count = popcnt64(sieve_[start_idx] & (m1 & m2));
  else
  {
    bit_count = popcnt64(sieve_[start_idx] & m1);
    bit_count += popcnt(&sieve_[start_idx + 1], stop_idx - (start_idx + 1));
    bit_count += popcnt64(sieve_[stop_idx] & m2);
  }

  return bit_count;
}

} // namespace
