///
/// @file  Sieve_pre_sieve.hpp
/// @brief Pre-sieve the (primes and) multiples of primes ≤ 71.
///        There are 7 static pre_sieved arrays in Sieve_arrays.hpp
///        from which the primes and multiples of primes have
///        been removed upfront. Each pre_sieved array corresponds
///        to exactly two sieving primes:
///
///        pre_sieved_arrays[0] = { 17, 19 }
///        pre_sieved_arrays[1] = { 23, 29 }
///        pre_sieved_arrays[2] = { 31, 37 }
///        pre_sieved_arrays[3] = { 41, 43 }
///        pre_sieved_arrays[4] = { 47, 53 }
///        pre_sieved_arrays[5] = { 59, 61 }
///        pre_sieved_arrays[6] = { 67, 71 }
///
///        Pre-sieving consists of bitwise AND'ing the values of
///        those 7 pre_sieved arrays and storing the result into the
///        the sieve array. Pre-sieving speeds up the S2_hard and
///        D algorithms by up to 5%.
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef SIEVE_PRE_SIEVE_HPP
#define SIEVE_PRE_SIEVE_HPP

#include <Sieve.hpp>
#include "Sieve_arrays.hpp"
#include <Vector.hpp>

#include <stdint.h>
#include <algorithm>
#include <cstring>

namespace {

/// Removes the (primes and) multiples of
/// primes ≤ 13 from the sieve array.
///
void pre_sieve1(primecount::Vector<uint8_t>& sieve, uint64_t low)
{
  uint64_t prime_product = pre_sieved_13.size() * 30;
  uint64_t i = (low % prime_product) / 30;
  uint64_t bytes_to_copy = pre_sieved_13.size() - i;
  const uint8_t* p = pre_sieved_13.begin();

  if (sieve.size() <= bytes_to_copy)
    std::copy_n(&p[i], sieve.size(), &sieve[0]);
  else
  {
    // Copy the last remaining bytes to the
    // beginning of the sieve array.
    std::copy_n(&p[i], bytes_to_copy, &sieve[0]);

    // Restart copying from the beginning
    for (i = bytes_to_copy; 
         i + pre_sieved_13.size() < sieve.size();
         i += pre_sieved_13.size())
      std::copy_n(p, pre_sieved_13.size(), &sieve[i]);

    // Copy the last remaining bytes
    std::copy_n(p, sieve.size() - i, &sieve[i]);
  }
}

void pre_sieve2(uint8_t* __restrict sieve,
                const uint8_t* __restrict pre_sieved,
                std::size_t bytes)
{
  constexpr std::size_t word_size = sizeof(unsigned long long);
  std::size_t limit = bytes - bytes % word_size;

  // Bitwise AND multiple bytes per iter.
  // std::memcpy is required to avoid unaligned
  // memory accesses. All major C++ compilers
  // will optimize away std::memcpy. 
  for (std::size_t i = 0; i < limit; i += word_size)
  {
    unsigned long long a, b;
    std::memcpy(&a, &sieve[i], word_size);
    std::memcpy(&b, &pre_sieved[i], word_size);
    unsigned long long result = a & b;
    std::memcpy(&sieve[i], &result, word_size);
  }

  // Bitwise AND the remaining bytes
  for (std::size_t i = limit; i < bytes; i++)
    sieve[i] &= pre_sieved[i];
}

} // namespace

namespace primecount {

/// Removes the (primes and) multiples of
/// primes ≤ 71 from the sieve array.
///
uint64_t Sieve::pre_sieve(uint64_t c, uint64_t low)
{
  // PrimePi(5) = 3
  uint64_t primePi = 3;

  if (c < 6)
    std::fill_n(sieve_.data(), sieve_.size(), 0xff);
  else
  {
    pre_sieve1(sieve_, low);
    // PrimePi(13) = 6
    primePi = 6;

    for (const auto& pre_sieved : pre_sieved_arrays)
    {
      // Each pre_sieved_arrays[i] contains
      // exactly 2 sieving primes > 13.
      if (c >= primePi + 2)
        primePi += 2;
      else
        return primePi;

      const uint8_t* p = pre_sieved.begin();
      uint64_t prime_product = pre_sieved.size() * 30;
      uint64_t pos = (low % prime_product) / 30;
      uint64_t offset = 0;

      while (offset < sieve_.size())
      {
        uint64_t bytes_to_copy = sieve_.size() - offset;
        uint64_t bytes_to_copy2 = uint64_t(pre_sieved.size() - pos);
        bytes_to_copy = std::min(bytes_to_copy, bytes_to_copy2);

        pre_sieve2(&sieve_[offset], &p[pos], bytes_to_copy);

        offset += bytes_to_copy;
        pos += bytes_to_copy;
        pos *= pos < pre_sieved.size();
      }
    }
  }

  return primePi;
}

} // namespace

#endif
