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
///         1) 0      if moebius(n) = 0
///         2) lpf    if !is_prime(n) && moebius(n) = -1
///         3) lpf-1  if !is_prime(n) && moebius(n) = 1
///         4) n      if  is_prime(n) && n < T_MAX
///         5) T_MAX  if  is_prime(n) && n > T_MAX
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef FACTORTABLE_HPP
#define FACTORTABLE_HPP

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <pmath.hpp>
#include <int128.hpp>

#include <algorithm>
#include <cassert>
#include <limits>
#include <stdint.h>
#include <vector>
#include <iostream>

#ifdef _OPENMP
  #include <omp.h>
#endif

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
  static int64_t get_index(uint64_t number)
  {
    assert(number > 0);
    uint64_t q = number / 210;
    uint64_t r = number % 210;
    return 48 * q + indexes_[r];
  }

  static int64_t get_number(uint64_t index)
  {
    uint64_t q = index / 48;
    uint64_t r = index % 48;
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
  /// Factor numbers <= y
  FactorTable(int64_t y, int threads)
  {
    if (y > max())
      throw primecount_error("y must be <= FactorTable::max().");
    y = std::max<int64_t>(8, y);
    T T_MAX = std::numeric_limits<T>::max();
    init_factors(y, T_MAX, threads);
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
  void init_factors(int64_t y, T T_MAX, int threads)
  {
    int64_t sqrty = isqrt(y);
    factors_.resize(get_index(y) + 1, T_MAX);
    factors_[0] = T_MAX - 1;

    int64_t thread_threshold = ipow(10, 7);
    threads = validate_threads(threads, y, thread_threshold);
    int64_t thread_interval = ceil_div(y, threads);

    #pragma omp parallel for num_threads(threads)
    for (int t = 0; t < threads; t++)
    {
      int64_t low = thread_interval * t;
      int64_t high = std::min(low + thread_interval, y);
      int64_t wheel_prime = get_number(1);
      primesieve::iterator iter(wheel_prime - 1, high);
      int64_t prime, i;

      while ((prime = iter.next_prime()) <= high)
      {
        // case 4), store prime
        if (prime > low && 
            prime < T_MAX)
          factors_[get_index(prime)] = (T) prime;

        int64_t multiple = next_multiple(low, prime, 1, &i);

        // no more work todo
        if (prime > T_MAX && multiple > high)
          break;

        for (; multiple <= high; multiple = prime * get_number(i++))
        {
          int64_t mi = get_index(multiple);
          // case 5), prime is the smallest factor of multiple
          if (factors_[mi] == T_MAX)
            factors_[mi] = (T) prime;
          // case 2) & 3), the least significant bit indicates
          // whether multiple has an even (0) or odd (1)
          // number of prime factors
          else if (factors_[mi] != 0)
            factors_[mi] ^= 1;
        }

        if (prime <= sqrty)
        {
          int64_t square = prime * prime;
          multiple = next_multiple(low, square, 0, &i);

          // case 1), set 0 if moebius(n) = 0
          for (; multiple <= high; multiple = square * get_number(i++))
            factors_[get_index(multiple)] = 0;
        }
      }
    }
  }

  /// Find the first multiple > low of prime which
  /// is not divisible by any prime <= 7.
  ///
  static int64_t next_multiple(int64_t low,
                               int64_t prime,
                               int64_t min_index,
                               int64_t* index)
  {
    int64_t quotient = (low / prime) + 1;
    int64_t i = std::max(min_index, get_index(quotient));
    int64_t multiple = prime * get_number(i++);

    for (; multiple <= low; i++)
      multiple = prime * get_number(i);
 
    *index = i;
    return multiple;
  }

  std::vector<T> factors_;
};

} // namespace

#endif
