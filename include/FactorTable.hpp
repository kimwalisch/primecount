///
/// @file  FactorTable.hpp
/// @brief The FactorTable class combines the lpf[n] (least prime
///        factor) and mu[n] (Möbius function) lookup tables into
///        a single factor_[n] table which furthermore only
///        contains entries for numbers which are not divisible by
///        2, 3, 5 and 7. The factor_[n] lookup table uses 17.5
///        times less memory than the lpf[n] & mu[n] lookup tables!
///        factor_[n] requires only 2 bytes per entry for 32-bit
///        numbers and 4 bytes per entry for 64-bit numbers.
///
///        The factor table concept has first been devised and
///        implemented by Christian Bau in 2003.
///
///        What we store in the factor_[n] lookup table:
///
///        1) INT_MAX      if n = 1
///        2) INT_MAX      if n is a prime
///        3) 0            if moebius(n) = 0
///        4) lpf - 1      if moebius(n) = 1
///        5) lpf          if moebius(n) = -1
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef FACTORTABLE_HPP
#define FACTORTABLE_HPP

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <imath.hpp>
#include <int128_t.hpp>

#include <algorithm>
#include <cassert>
#include <limits>
#include <stdint.h>
#include <vector>

namespace primecount {

/// AbstractFactorTable contains static lookup tables
/// and is used to convert:
/// 1) A number into a FactorTable index
/// 2) A FactorTable index into a number
///
class AbstractFactorTable
{
public:
  static void to_index(int64_t* number)
  {
    assert(*number > 0);
    *number = get_index(*number);
  }

  static int64_t get_index(uint64_t number)
  {
    assert(number > 0);
    uint64_t q = number / 2310;
    uint64_t r = number % 2310;
    return 480 * q + indexes_[r];
  }

  static int64_t get_number(uint64_t index)
  {
    uint64_t q = index / 480;
    uint64_t r = index % 480;
    return 2310 * q + numbers_[r];
  }

private:
  static const uint16_t numbers_[480];
  static const int16_t indexes_[2310];
};

template <typename T>
class FactorTable : public AbstractFactorTable
{
public:
  /// Factor numbers <= y
  FactorTable(int64_t y, int threads)
  {
    if (y > max())
      throw primecount_error("y must be <= FactorTable::max()");

    y = std::max<int64_t>(8, y);
    T T_MAX = std::numeric_limits<T>::max();
    factor_.resize(get_index(y) + 1, T_MAX);

    int64_t sqrty = isqrt(y);
    int64_t thread_threshold = ipow(10, 7);
    threads = ideal_num_threads(threads, y, thread_threshold);
    int64_t thread_distance = ceil_div(y, threads);

    #pragma omp parallel for num_threads(threads)
    for (int t = 0; t < threads; t++)
    {
      int64_t low = 1;
      low += thread_distance * t;
      int64_t high = std::min(low + thread_distance, y);
      primesieve::iterator it(get_number(1) - 1);

      while (true)
      {
        int64_t i = 1;
        int64_t prime = it.next_prime();
        int64_t multiple = next_multiple(prime, low, &i);
        int64_t min_m = prime * get_number(1);

        if (min_m > high)
          break;

        for (; multiple <= high; multiple = prime * get_number(i++))
        {
          int64_t mi = get_index(multiple);
          // prime is smallest factor of multiple
          if (factor_[mi] == T_MAX)
            factor_[mi] = (T) prime;
          // the least significant bit indicates
          // whether multiple has an even (0) or odd (1)
          // number of prime factors
          else if (factor_[mi] != 0)
            factor_[mi] ^= 1;
        }

        if (prime <= sqrty)
        {
          int64_t j = 0;
          int64_t square = prime * prime;
          multiple = next_multiple(square, low, &j);

          // moebius(n) = 0
          for (; multiple <= high; multiple = square * get_number(j++))
            factor_[get_index(multiple)] = 0;
        }
      }
    }
  }

  /// Get the least prime factor (lpf) of the number
  /// n = get_number(index). The return value is different
  /// from the least prime factor in some situations
  /// but this does not affect our calculations.
  ///
  /// 1) INT_MAX      if n = 1
  /// 2) INT_MAX      if n is a prime
  /// 3) 0            if moebius(n) = 0
  /// 4) lpf - 1      if moebius(n) = 1
  /// 5) lpf          if moebius(n) = -1
  ///
  int64_t lpf(int64_t index) const
  {
    return factor_[index];
  }

  /// Get the Möbius function value of the number
  /// n = get_number(index). For performance reasons
  /// mu(index) == 0 is not supported.
  ///
  int64_t mu(int64_t index) const
  {
    assert(factor_[index] != 0);
    return (factor_[index] & 1) ? -1 : 1;
  }

  static maxint_t max()
  {
    maxint_t T_MAX = std::numeric_limits<T>::max();
    return ipow(T_MAX - 1, 2) - 1;
  }

private:

  /// Find the first multiple (of prime) > low which
  /// is not divisible by any prime <= 7
  ///
  static int64_t next_multiple(int64_t prime,
                               int64_t low,
                               int64_t* index)
  {
    int64_t quotient = ceil_div(low, prime);
    int64_t i = std::max(*index, get_index(quotient));
    int64_t multiple = 0;

    for (; multiple <= low; i++)
      multiple = prime * get_number(i);

    *index = i;
    return multiple;
  }

  std::vector<T> factor_;
};

} // namespace

#endif
