///
/// @file  bit_sieve.hpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef BIT_SIEVE_HPP
#define BIT_SIEVE_HPP

#include <primecount-internal.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <vector>
#include <stdint.h>

namespace primecount {

/// The bit_sieve data structure uses bit packing to save memory.
/// bit_sieve uses 1 byte of memory for 16 numbers, each bit
/// corresponds to one odd integer.
///
class bit_sieve
{
public:
  bit_sieve(std::size_t size)
    : bits_((size + 16 - 1) / 16 + sizeof(uint64_t)),
      size_(size),
      low_(0)
  { }

  /// Count the number of 1 bits inside [start, stop]
  int64_t count(int64_t start, int64_t stop) const
  {
    assert(stop < size_);
    const uint64_t* bits = reinterpret_cast<const uint64_t*>(bits_.data());
    return popcount(bits, start, stop, low_);
  }

  /// Set all bits to 1
  void fill()
  {
    std::fill(bits_.begin(), bits_.end(), 0xff);

    // As bits correspond to odd integers we must
    // recalculate the bitmasks if low changes in
    // order to ensure bit_sieve[even] = false
    //
    for (int i = 0; i < 16; i += 2)
    {
      int mask0 = 1 << (i / 2);
      int mask1 = 0;

      if (low_ % 2 == 0)
        std::swap(mask0, mask1);

      is_bit[i + 0] = mask0;
      is_bit[i + 1] = mask1;

      unset_bit[i + 0] = ~mask0;
      unset_bit[i + 1] = ~mask1;
    }
  }

  bool operator[](uint64_t pos) const
  {
    assert(pos < size_);
    return (bits_[pos >> 4] & is_bit[pos & 15]) != 0;
  }

  std::size_t size() const
  {
    return size_;
  }

  void set_low(int64_t low)
  {
    low_ = low;
  }

  void unset(uint64_t pos)
  {
    assert(pos < size_);
    bits_[pos >> 4] &= unset_bit[pos & 15];
  }
private:
  int is_bit[16];
  int unset_bit[16];
  std::vector<unsigned char> bits_;
  std::size_t size_;
  int64_t low_;
};

} // namespace

#endif
