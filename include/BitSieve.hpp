///
/// @file  BitSieve.hpp
/// @brief The BitSieve class is a bit array for prime sieving
///        that packs 64 numbers into 8 bytes i.e. each bit
///        corresponds to one integer.
///
/// Copyright (C) 2015 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef BITSIEVE_HPP
#define BITSIEVE_HPP

#include <cassert>
#include <cstddef>
#include <vector>
#include <stdint.h>

namespace primecount {

class BitSieve
{
public:
  BitSieve(std::size_t size);

  /// Set all bits to 1, except bits corresponding to even numbers.
  /// @warning You must reset 2 when sieving primes.
  ///
  void fill(uint64_t low, uint64_t /* unused */);

  /// Sieve out the multiples of the first c primes.
  void pre_sieve(uint64_t c, uint64_t low);

  /// Count the number of 1 bits inside the interval [start, stop]
  uint64_t count(uint64_t start, uint64_t stop) const;

  /// Count the number of 1 bits inside [0, stop]
  uint64_t count(uint64_t stop) const
  {
    return count(0, stop);
  }

  /// Count the number of 1 bits inside the interval [start, stop]
  uint64_t count(uint64_t start,
                 uint64_t stop,
                 uint64_t low,
                 uint64_t high,
                 uint64_t count_0_start,
                 uint64_t count_low_high) const
  {
    if (start > stop)
      return 0;

    if (stop - start < high - low - stop)
      return count(start, stop);
    else
      // optimization, same as count(start, stop)
      return count_low_high - count_0_start - count(stop + 1, (high - 1) - low);
  }

  bool operator[](uint64_t pos) const
  {
    assert(pos < size_);
    uint64_t bit = ((uint64_t) 1) << (pos & 63);
    return (bits_[pos >> 6] & bit) != 0;
  }

  void set(uint64_t pos)
  {
    assert(pos < size_);
    uint64_t bit = ((uint64_t) 1) << (pos & 63);
    bits_[pos >> 6] |= bit;
  }

  void unset(uint64_t pos)
  {
    assert(pos < size_);
    bits_[pos >> 6] &= unset_bit_[pos & 63];
  }

  std::size_t size() const
  {
    return size_;
  }
private:
  static const uint64_t unset_bit_[64];
  std::vector<uint64_t> bits_;
  std::size_t size_;
};

} // namespace

#endif
