///
/// @file  FactorTable.hpp
/// @brief The FactorTable class is used to save memory. It combines
///        the lpf[n] (least prime factor) and mu[n] (Möbius function)
///        lookup tables into a single factors_[n] table which
///        furthermore only contains entries for numbers which are
///        not divisible by 2, 3, 5 and 7.
///        The factor table concept has first been devised and
///        implemented by Christian Bau in 2003.
///
///        FactorTable.lpf(index) is equal to (n = get_number(index)):
///
///         * 0      if moebius(n) = 0
///         * lpf    if !is_prime(n) && moebius(n) = -1
///         * lpf-1  if !is_prime(n) && moebius(n) = 1
///         * n      if  is_prime(n) && n < T_MAX
///         * T_MAX  if  is_prime(n) && n > T_MAX
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef FACTORTABLE_HPP
#define FACTORTABLE_HPP

#include <primecount.hpp>
#include <pmath.hpp>
#include <int128.hpp>

#include <algorithm>
#include <cassert>
#include <limits>
#include <stdint.h>
#include <vector>

namespace primecount {

/// AbstractFactorTable contains static FactorTable
/// data and it used to convert:
/// 1) A number into a FactorTable index.
/// 2) A FactorTable index into a number.
///
class AbstractFactorTable
{
protected:
  virtual ~AbstractFactorTable() { }

public:
  /// @pre number > 0
  static void to_index(int64_t* number)
  {
    assert(*number > 0);
    *number = get_index(*number);
  }

  /// @pre number > 0
  static int64_t get_index(int64_t number)
  {
    assert(number > 0);
    int64_t q = ((uint64_t) number) / 210;
    int64_t r = number - q * 210;

    return 48 * q + indexes_[r];
  }

  static int64_t get_number(int64_t index)
  {
    // unsigned is faster here
    int64_t q = ((uint64_t) index) / 48;
    int64_t r = index - q * 48;

    return 210 * q + numbers_[r];
  }

private:
  static const uint8_t numbers_[48];
  static const  int8_t indexes_[210];
};

template <typename T>
class FactorTable : public AbstractFactorTable
{
public:
  /// @param y  factor numbers <= y
  FactorTable(int64_t y)
  {
    if (y > max())
      throw primecount_error("y must be <= FactorTable::max().");
    y = std::max<int64_t>(8, y);
    T T_MAX = std::numeric_limits<T>::max();
    init_factors(y, T_MAX);
  }

  static maxint_t max()
  {
    maxint_t T_MAX = std::numeric_limits<T>::max();
    return ipow(T_MAX - 1, 2) - 1;
  }

  /// Get the least prime factor (lpf) of the number get_number(index).
  /// The result is different from lpf in some situations:
  /// 1) lpf(index) returns T_MAX if get_number(index) is a prime > T_MAX.
  /// 2) lpf(index) returns lpf minus one if mu(index) == 1.
  /// 3) lpf(index) returns 0 if get_number(index) has a squared prime factor.
  ///
  int64_t lpf(int64_t index) const
  {
    return factors_[index];
  }

  /// Get the Möbius function value of the number get_number(index).
  /// For performance reasons mu(index) == 0 is not supported.
  /// @pre mu(index) != 0 (i.e. lpf(index) != 0)
  ///
  int64_t mu(int64_t index) const
  {
    assert(lpf(index) != 0);
    return (factors_[index] & 1) ? -1 : 1;
  }

private:
  void init_factors(int64_t y, T T_MAX)
  {
    int64_t sqrty = isqrt(y);
    factors_.resize(get_index(y) + 1, T_MAX);
    factors_[0] = T_MAX - 1;

    for (std::size_t i = 1; i < factors_.size(); i++)
    {
      if (factors_[i] == T_MAX)
      {
        int64_t prime = get_number(i);
        int64_t multiple = prime * get_number(1);

        if (prime < T_MAX)
          factors_[i] = (T) prime;
        else if (multiple > y)
          break;

        for (int64_t j = 2; multiple <= y; multiple = prime * get_number(j++))
        {
          int64_t index = get_index(multiple);
          // prime is the smallest factor of multiple
          if (factors_[index] == T_MAX)
            factors_[index] = (T) prime;
          // the least significant bit indicates whether multiple has
          // an even (0) or odd (1) number of prime factors
          else if (factors_[index] != 0)
            factors_[index] ^= 1;
        }

        // Möbius function is 0 if n has a squared prime factor
        if (prime <= sqrty)
        {
          multiple = prime * prime;
          for (int64_t j = 1; multiple <= y; multiple = prime * prime * get_number(j++))
            factors_[get_index(multiple)] = 0;
        }
      }
    }
  }

  std::vector<T> factors_;
};

} // namespace

#endif
