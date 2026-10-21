///
/// @file  FactorTable.hpp
/// @brief The FactorTable class combines the lpf[n] (least prime
///        factor) and mu[n] (Möbius function) lookup tables into a
///        single factor[n] table which furthermore only contains
///        entries for numbers which are not divisible by 2, 3, 5, 7
///        and 11. The factor[n] lookup table uses up to 12.03
///        times less memory than the lpf[n] & mu[n] lookup tables!
///        factor[n] uses only 2 or 4 bytes per entry.
///
///        The factor table concept was devised and implemented by
///        Christian Bau in 2003. Note that Tomás Oliveira e Silva
///        also suggests combining the mu[n] and lpf[n] lookup tables
///        in his paper. However Christian Bau's FactorTable data
///        structure uses only half as much memory and is also
///        slightly more efficient (uses fewer instructions) than the
///        data structure proposed by Tomás Oliveira e Silva.
///
///        factor[n] stores the index pi(lpf) of the least prime
///        factor and the Möbius sign. As pi(lpf) is far smaller
///        than lpf this allows 2 byte entries up to a much larger n:
///
///        1) INT_MAX - 1    if n = 1
///        2) INT_MAX        if n is a prime
///        3) 0              if moebius(n) = 0
///        4) 2*pi(lpf)      if moebius(n) = 1
///        5) 2*pi(lpf) + 1  if moebius(n) = -1
///
///        factor[1] = (INT_MAX - 1) because 1 contributes to the
///        sum of the ordinary leaves S1(x, a) in the
///        Lagarias-Miller-Odlyzko and Deleglise-Rivat algorithms.
///        The values above allow to replace the 1st if statement
///        below, used in the S2_hard (Deleglise-Rivat) and
///        D (Gourdon) formulas, by the 2nd new if statement which
///        is obviously faster.
///
///        * Old: if (mu[n] != 0 && lpf[n] > prime)
///        * New: if (factor[n] > 2*pi(prime) + 1)
///
///        In-depth description of the factor table data structure:
///        https://github.com/kimwalisch/primecount/blob/master/doc/Hard-Special-Leaves-SIMD-Filtering.pdf
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef FACTORTABLE_HPP
#define FACTORTABLE_HPP

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <BaseFactorTable.hpp>
#include <primesieve.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <macros.hpp>
#include <Vector.hpp>

#include <algorithm>
#include <stdint.h>

namespace {

using namespace primecount;

template <typename T>
class FactorTable : public BaseFactorTable
{
public:
  static_assert(sizeof(T) == sizeof(uint16_t) ||
                sizeof(T) == sizeof(uint32_t),
                "FactorTable: T must be uint16_t or uint32_t!");

  /// Factor numbers <= y
  FactorTable(int64_t y, int threads)
  {
    if_unlikely(y > max())
      throw primecount_error("y must be <= FactorTable::max()");

    y = std::max<int64_t>(1, y);
    T T_MAX = pstd::numeric_limits<T>::max();
    factor_.resize(to_index(y) + 1);

    // mu(1) = 1.
    // 1 has zero prime factors, hence 1 has an even
    // number of prime factors. We use the least
    // significant bit to indicate whether the number
    // has an even or odd number of prime factors.
    factor_[0] = T_MAX ^ 1;

    int64_t sqrty = isqrt(y);
    int64_t thread_threshold = (int64_t) 5e6;
    threads = ideal_num_threads(y, threads, thread_threshold);
    int64_t thread_dist = ceil_div(y, threads);
    thread_dist += coprime_indexes_.size() - thread_dist % coprime_indexes_.size();

    #pragma omp parallel for num_threads(threads) schedule(static, 1)
    for (int t = 0; t < threads; t++)
    {
      // Thread processes interval [low, high]
      int64_t low = thread_dist * t;
      int64_t high = low + thread_dist;
      int64_t min_m = first_coprime() * first_coprime();
      low = std::max(first_coprime(), low + 1);
      high = std::min(high, y);

      if (low <= high &&
          min_m <= high)
      {
        // Default initialize memory to all bits set
        int64_t low_idx = to_index(low);
        int64_t size = (to_index(high) + 1) - low_idx;
        std::fill_n(&factor_[low_idx], size, T_MAX);

        int64_t start = first_coprime();
        int64_t stop = high / first_coprime();
        primesieve::iterator it(start, stop);

        // PrimePi(13) = 6
        ASSERT(first_coprime() == 13);
        int64_t prime_index = 6;

        for (; true; prime_index++)
        {
          int64_t i = 1;
          int64_t prime = it.next_prime();
          int64_t multiple = next_multiple(prime, low, &i);
          min_m = prime * first_coprime();

          if (min_m > high)
            break;

          for (; multiple <= high; multiple = prime * to_number(i++))
          {
            int64_t mi = to_index(multiple);
            // prime is smallest factor of multiple
            if (factor_[mi] == T_MAX)
            {
              ASSERT(encode(prime_index) <= T_MAX - 2);
              factor_[mi] = (T) encode(prime_index);
            }
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
            for (; multiple <= high; multiple = square * to_number(j++))
              factor_[to_index(multiple)] = 0;
          }
        }
      }
    }
  }

  /// Returns the factor table entry for the number
  /// n = to_number(index).
  ///
  /// Return value:
  ///
  /// 1) INT_MAX - 1    if n = 1
  /// 2) INT_MAX        if n is a prime
  /// 3) 0              if moebius(n) = 0
  /// 4) 2*pi(lpf)      if moebius(n) = 1
  /// 5) 2*pi(lpf) + 1  if moebius(n) = -1
  ///
  int64_t operator[](int64_t index) const
  {
    return factor_[index];
  }

  const T* data() const
  {
    return factor_.data();
  }

  /// Encode prime index for use in the factor table
  static int64_t encode(int64_t prime_index)
  {
    return prime_index * 2 + 1;
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

  static constexpr int64_t max()
  {
    // The least prime factor is <= sqrt(y) and with 16-bit
    // entries p(32767) = 386083 is the first prime whose
    // index cannot be encoded, hence y < 386083^2.
    return sizeof(T) == sizeof(uint16_t)
      ? 386083ll * 386083 - 1
      : pstd::numeric_limits<int64_t>::max();
  }

private:
  Vector<T> factor_;
};

} // namespace

#endif
