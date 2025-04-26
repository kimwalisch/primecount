///
/// @file  Sieve_pre_sieve.hpp
/// @brief Pre-sieve the primes and multiples of primes ≤ 37.
///        There are 3 static pre_sieved arrays in Sieve_arrays.hpp
///        from which the primes and multiples of primes have
///        been removed upfront. Each pre_sieved array corresponds
///        to exactly three sieving primes:
///
///        pre_sieved_arrays[0] = {  7, 11, 13 }
///        pre_sieved_arrays[1] = { 17, 19, 23 }
///        pre_sieved_arrays[2] = { 29, 31, 37 }
///
///        Pre-sieving consists of bitwise AND'ing the values of
///        those 3 pre_sieved arrays and storing the result into the
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

/// Initializes the sieve array 
void pre_sieve1(primecount::Vector<uint8_t>& sieve, uint64_t low)
{
  const auto& pre_sieved = pre_sieved_arrays[0];
  const uint8_t* ptr = &*pre_sieved.begin();
  uint64_t prime_product = pre_sieved.size() * 30;
  uint64_t i = (low % prime_product) / 30;
  uint64_t bytes_to_copy = pre_sieved.size() - i;

  if (sieve.size() <= bytes_to_copy)
    std::copy_n(&ptr[i], sieve.size(), &sieve[0]);
  else
  {
    // Copy the last remaining bytes to the
    // beginning of the sieve array.
    std::copy_n(&ptr[i], bytes_to_copy, &sieve[0]);

    // Restart copying from the beginning
    for (i = bytes_to_copy; 
         i + pre_sieved.size() < sieve.size();
         i += pre_sieved.size())
      std::copy_n(ptr, pre_sieved.size(), &sieve[i]);

    // Copy the last remaining bytes
    std::copy_n(ptr, sieve.size() - i, &sieve[i]);
  }
}

void pre_sieve2(uint8_t* __restrict sieve,
                const uint8_t* __restrict pre_sieved,
                std::size_t bytes)
{
  constexpr std::size_t word_size = sizeof(unsigned long long);
  std::size_t limit = bytes - bytes % word_size;

  // Bitwise AND multiple bytes per iter
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

uint64_t Sieve::pre_sieve(uint64_t c, uint64_t low)
{
  // PrimePi(x) of max prime used
  // in pre_sieved_arrays[i].
  uint64_t c2 = 3;

  if (c < c2 + 3)
    std::fill_n(sieve_.data(), sieve_.size(), 0xff);
  else
  {
    pre_sieve1(sieve_, low);

    for (std::size_t i = 1; i < pre_sieved_arrays.size(); i++)
    {
      // Each pre_sieved lookup table
      // contains exactly 3 primes.
      c2 += 3;
      if (c < c2 + 3)
        return c2;

      const auto& pre_sieved = pre_sieved_arrays[i];
      const uint8_t* ptr = &*pre_sieved.begin();
      uint64_t prime_product = pre_sieved.size() * 30;
      uint64_t pos = (low % prime_product) / 30;
      uint64_t offset = 0;

      while (offset < sieve_.size())
      {
        uint64_t bytes_to_copy = sieve_.size() - offset;
        uint64_t bytes_to_copy2 = uint64_t(pre_sieved.size() - pos);
        bytes_to_copy = std::min(bytes_to_copy, bytes_to_copy2);

        pre_sieve2(&sieve_[offset], &ptr[pos], bytes_to_copy);

        offset += bytes_to_copy;
        pos += bytes_to_copy;
        pos *= pos < pre_sieved.size();
      }
    }
  }

  return c2;
}

} // namespace

#endif
