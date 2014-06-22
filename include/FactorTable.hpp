///
/// @file  FactorTable.hpp
/// @brief The FactorTable class is used to save memory. It combines
///        the lpf[n] (least prime factor) and mu[n] (Möbius function)
///        lookup tables into a single factor_table[n] which
///        furthermore only contains entries for numbers which are
///        not divisible by 2, 3, 5 and 7.
///        The factor table concept has first been devised and
///        implemented by Christian Bau in 2003.
///
///        FactorTable.lpf(index) is equal to (n = get_number(index)):
///
///         * 0      if moebius(n) = 0
///         * lpf    if !is_prime(n) && moebius(n) = 1
///         * lpf-1  if !is_prime(n) && moebius(n) = -1
///         * n      if  is_prime(n) && n < 65535
///         * 65535  if  is_prime(n) && n > 65535
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef FACTORTABLE_HPP
#define FACTORTABLE_HPP

#include <cassert>
#include <stdint.h>
#include <vector>

class FactorTable
{
public:
  FactorTable(int64_t max, int64_t sqrty);
  void init();

  /// Get the FactorTable index corressponding to the number.
  /// @pre get_number(index) > 7
  ///
  int64_t get_index(int64_t number) const
  {
    assert(get_number(index) > 7);
    return 48 * (number / 210) + indexes_[number % 210];
  }

  /// Get the number corresponding to the FactorTable index.
  /// @pre get_number(index) > 7
  ///
  int64_t get_number(int64_t index) const
  {
    assert(get_number(index) > 7);
    return 210 * (index / 48) + numbers_[index % 48];
  }

  /// Get the least prime factor of the number get_number(index).
  /// Note that if get_number(index) is a prime > 65535 than
  /// lpf(index) returns 65535 and if mu(index) == -1 then lpf(index)
  /// returns the least prime factor minus one.
  /// @pre get_number(index) > 7
  ///
  int64_t lpf(int64_t index) const
  {
    assert(get_number(index) > 7);
    return factor_table_[index];
  }

  /// Get the Möbius function value of the number get_number(index).
  /// @pre lpf(index) != 0
  ///
  int64_t mu(int64_t index) const
  {
    assert(lpf(index) != 0);
    return (factor_table_[index] & 1) ? 1 : -1;
  }
private:
  static const  int8_t indexes_[210];
  static const uint8_t numbers_[48];

  std::vector<uint16_t> factor_table_;
  int64_t max_;
};

#endif
