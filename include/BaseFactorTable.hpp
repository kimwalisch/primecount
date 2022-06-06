///
/// @file  BaseFactorTable.hpp
///        Static lookup tables and functions used by the
///        FactorTable and FactorTableD classes.
///        See FactorTable.hpp for more information.
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef BASEFACTORTABLE_HPP
#define BASEFACTORTABLE_HPP

#include <imath.hpp>
#include <macros.hpp>
#include <pod_vector.hpp>

#include <algorithm>
#include <stdint.h>

namespace primecount {

/// BaseFactorTable contains static lookup tables
/// and is used to convert:
/// 1) A number into a FactorTable index
/// 2) A FactorTable index into a number
///
class BaseFactorTable
{
public:
  static int64_t to_index(uint64_t number)
  {
    ASSERT(number > 0);
    uint64_t q = number / 2310;
    uint64_t r = number % 2310;
    return 480 * q + coprime_indexes_[r];
  }

  static int64_t to_number(uint64_t index)
  {
    uint64_t q = index / 480;
    uint64_t r = index % 480;
    return 2310 * q + coprime_[r];
  }

  /// Returns the 1st number > 1 that is not divisible
  /// by 2, 3, 5, 7 and 11. Hence 13 is returned.
  ///
  static int64_t first_coprime()
  {
    return to_number(1);
  }

protected:
  /// Find the first multiple (of prime) >= low which
  /// is not divisible by any prime <= 11.
  ///
  static int64_t next_multiple(int64_t prime,
                               int64_t low,
                               int64_t* index)
  {
    int64_t quotient = ceil_div(low, prime);
    int64_t i = std::max(*index, to_index(quotient));
    int64_t multiple = 0;

    for (; multiple < low; i++)
      multiple = prime * to_number(i);

    *index = i;
    return multiple;
  }

  static const pod_array<uint16_t, 480> coprime_;
  static const pod_array<int16_t, 2310> coprime_indexes_;
};

} // namespace 

#endif
