///
/// @file  FactorTableD.hpp
/// @brief The FactorTableD class combines the lpf[n] (least prime
///        factor), mpf[n] (max prime factor) and mu[n] (Möbius
///        function) lookup tables into a single factor[n] lookup
///        table which furthermore only contains entries for numbers
///        which are not divisible by 2, 3, 5, 7 and 11. The factor[n]
///        lookup table uses up to 28 times less memory than the
///        lpf[n], mpf[n] and mu[n] lookup tables! factor[n] uses only
///        2 bytes per entry for 32-bit numbers and 4 bytes per entry
///        for 64-bit numbers.
///
///        The factor table concept was devised and implemented by
///        Christian Bau in 2003. Note that Tomás Oliveira e Silva
///        also suggests combining the mu[n] and lpf[n] lookup tables
///        in his paper. However Christian Bau's FactorTable data
///        structure uses only half as much memory and is also
///        slightly more efficient (uses fewer instructions) than the
///        data structure proposed by Tomás Oliveira e Silva.
///
///        What we store in the factor[n] lookup table:
///
///        1) INT_MAX - 1  if n = 1
///        2) INT_MAX      if n is a prime
///        3) 0            if n has a prime factor > y
///        4) 0            if moebius(n) = 0
///        5) lpf - 1      if moebius(n) = 1
///        6) lpf          if moebius(n) = -1
///
///        factor[1] = (INT_MAX - 1) because 1 contributes to the
///        sum of the ordinary leaves S1(x, a) in the
///        Lagarias-Miller-Odlyzko and Deleglise-Rivat algorithms.
///        The values above allow to replace the 1st if statement
///        below used in the D(x, y) formula by the 2nd new if
///        statement which is obviously faster.
///
///        * Old: if (mu[n] != 0 && lpf[n] > prime && mpf[n] <= y)
///        * New: if (prime < factor[n])
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef FACTORTABLED_HPP
#define FACTORTABLED_HPP

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <BaseFactorTable.hpp>
#include <primesieve.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <macros.hpp>
#include <pod_vector.hpp>

#include <algorithm>
#include <limits>
#include <stdint.h>

namespace {

using namespace primecount;

template <typename T>
class FactorTableD : public BaseFactorTable
{
public:
  /// Factor numbers <= z
  FactorTableD(int64_t y,
               int64_t z,
               int threads)
  {
    if (z > max())
      throw primecount_error("z must be <= FactorTable::max()");

    z = std::max<int64_t>(1, z);
    T T_MAX = std::numeric_limits<T>::max();
    factor_.resize(to_index(z) + 1);

    // mu(1) = 1.
    // 1 has zero prime factors, hence 1 has an even
    // number of prime factors. We use the least
    // significant bit to indicate whether the number
    // has an even or odd number of prime factors.
    factor_[0] = T_MAX ^ 1;

    int64_t sqrtz = isqrt(z);
    int64_t thread_threshold = (int64_t) 1e7;
    threads = ideal_num_threads(threads, z, thread_threshold);
    int64_t thread_distance = ceil_div(z, threads);
    thread_distance += coprime_indexes_.size() - thread_distance % coprime_indexes_.size();

    #pragma omp parallel for num_threads(threads)
    for (int t = 0; t < threads; t++)
    {
      // Thread processes interval [low, high]
      int64_t low = thread_distance * t;
      int64_t high = low + thread_distance;
      low = std::max(first_coprime(), low + 1);
      high = std::min(high, z);

      if (low <= high)
      {
        // Default initialize memory to all bits set
        int64_t low_idx = to_index(low);
        int64_t size = (to_index(high) + 1) - low_idx;
        std::fill_n(&factor_[low_idx], size, T_MAX);

        int64_t start = first_coprime() - 1;
        int64_t stop = high / first_coprime();
        int64_t min_m = first_coprime() * first_coprime();
        primesieve::iterator it(start, stop);

        if (min_m <= high)
        {
          while (true)
          {
            // Find multiples > prime
            int64_t i = 1;
            int64_t prime = it.next_prime();
            int64_t multiple = next_multiple(prime, low, &i);
            min_m = prime * first_coprime();

            if (min_m > high)
              break;

            for (; multiple <= high; multiple = prime * to_number(i++))
            {
              int64_t mi = to_index(multiple);
              // prime is the smallest factor of multiple
              if (factor_[mi] == T_MAX)
                factor_[mi] = (T) prime;
              // the least significant bit indicates
              // whether multiple has an even (0) or odd (1)
              // number of prime factors
              else if (factor_[mi] != 0)
                factor_[mi] ^= 1;
            }

            if (prime <= sqrtz)
            {
              int64_t j = 0;
              int64_t square = prime * prime;
              multiple = next_multiple(square, low, &j);

              // Sieve out numbers that are not square free
              // i.e. numbers for which moebius(n) = 0.
              for (; multiple <= high; multiple = square * to_number(j++))
                factor_[to_index(multiple)] = 0;
            }
          }
        }

        start = std::max(start, y);

        if (start < high)
        {
          it.skipto(start, high);

          // y < prime <= z
          while (true)
          {
            int64_t i = 0;
            int64_t prime = it.next_prime();
            int64_t next = next_multiple(prime, low, &i);

            if (prime > high)
              break;

            // Sieve out primes > y &&
            // Sieve out numbers with prime factors > y
            for (; next <= high; next = prime * to_number(i++))
              factor_[to_index(next)] = 0;
          }
        }
      }
    }
  }

  /// Returns true if n (with n = to_number(index)) is a
  /// hard special leaf in the D formula of Xavier
  /// Gourdon's prime counting algorithm.
  ///
  /// Return value:
  ///
  /// 1) INT_MAX - 1  if n = 1
  /// 2) INT_MAX      if n is a prime
  /// 3) 0            if n has a prime factor > y
  /// 4) 0            if moebius(n) = 0
  /// 5) lpf - 1      if moebius(n) = 1
  /// 6) lpf          if moebius(n) = -1
  ///
  int64_t is_leaf(int64_t index) const
  {
    return factor_[index];
  }

  /// Get the Möbius function value of the number
  /// n = to_number(index).
  ///
  /// https://en.wikipedia.org/wiki/Möbius_function
  /// mu(n) = 1 if n is a square-free integer with an even number of prime factors.
  /// mu(n) = −1 if n is a square-free integer with an odd number of prime factors.
  /// mu(n) = 0 if n has a squared prime factor.
  ///
  int64_t mu(int64_t index) const
  {
    // mu(n) = 0 is disabled by default for performance
    // reasons, we only enable it for testing.
    #if defined(ENABLE_MU_0_TESTING)
      if (factor_[index] == 0)
        return 0;
    #else
      ASSERT(factor_[index] != 0);
    #endif

    if (factor_[index] & 1)
      return -1;
    else
      return 1;
  }

  static maxint_t max()
  {
    maxint_t T_MAX = std::numeric_limits<T>::max();
    return ipow(T_MAX - 1, 2) - 1;
  }

private:
  pod_vector<T> factor_;
};

} // namespace

#endif
