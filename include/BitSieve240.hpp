///
/// @file  BitSieve240.hpp
/// @brief The BitSieve240 base class contains lookup tables that are
///        needed to implement a prime sieving algorithm where each
///        bit corresponds to an integer that is not divisible by 2,
///        3 and 5. The 8 bits of each byte correspond to the offsets
///        { 1, 7, 11, 13, 17, 19, 23, 29 }. Since the sieve array
///        uses the uint64_t data type, one sieve array element
///        (8 bytes) corresponds to an interval of size 30 * 8 = 240.
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef BITSIEVE240_HPP
#define BITSIEVE240_HPP

#include <pod_vector.hpp>
#include <stdint.h>

namespace primecount {

class BitSieve240
{
protected:
  static const pod_array<uint64_t, 6> pi_tiny_;
  static const pod_array<uint64_t, 240> set_bit_;
  static const pod_array<uint64_t, 240> unset_bit_;
  static const pod_array<uint64_t, 240> unset_larger_;
};

} // namespace

#endif
