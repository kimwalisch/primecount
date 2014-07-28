///
/// @file  S1.cpp
/// @brief Functions to calculate the contribution of the ordinary
///        leaves in the Lagarias-Miller-Odlyzko algorithm.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <FactorTable.hpp>
#include <PhiTiny.hpp>
#include <pmath.hpp>
#include <ptypes.hpp>

#include <stdint.h>
#include <vector>

using namespace std;
using namespace primecount;

namespace {
namespace S1 {

/// Run time: O(y) operations, O(y) space.
/// @pre is_phi_tiny(c) == true
///
template <typename T>
T S1(T x,
     int64_t y,
     int64_t c,
     vector<int32_t>& primes,
     vector<int32_t>& lpf,
     vector<int32_t>& mu)
{
  T S1_result = 0;

  for (int64_t n = 1; n <= y; n++)
    if (lpf[n] > primes[c])
      S1_result += mu[n] * phi_tiny(x / n, c);

  return S1_result;
}

/// This version uses 17 times less memory than the version above.
/// Run time: O(y) operations, O(y) space.
/// @pre is_phi_tiny(c) == true
///
template <typename T>
T S1(T x,
     int64_t y,
     int64_t c,
     vector<int32_t>& primes,
     FactorTable& factors)
{
  // FactorTable contains only numbers coprime to 2, 3, 5 and 7
  if (primes[c] <= 7)
  {
    vector<int32_t> mu = generate_moebius(y);
    vector<int32_t> lpf = generate_least_prime_factors(y);
    return S1(x, y, c, primes, lpf, mu);
  }

  T S1_result = 0;
  int64_t limit = factors.get_index(y);

  for (int64_t i = factors.get_index(1); i <= limit; i++)
    if (factors.lpf(i) > primes[c])
      S1_result += factors.mu(i) * phi_tiny(x / factors.get_number(i), c);

  return S1_result;
}

} // namespace
} // namespace S1

namespace primecount {

int64_t S1(int64_t x,
           int64_t y,
           int64_t c,
           vector<int32_t>& primes,
           vector<int32_t>& lpf,
           vector<int32_t>& mu)
{
  return S1::S1(x, y, c, primes, lpf, mu);
}

int64_t S1(int64_t x,
           int64_t y,
           int64_t c,
           vector<int32_t>& primes,
           FactorTable& factors)
{
  return S1::S1(x, y, c, primes, factors);
}

#ifdef HAVE_INT128_T

int128_t S1(int128_t x,
            int64_t y,
            int64_t c,
            vector<int32_t>& primes,
            vector<int32_t>& lpf,
            vector<int32_t>& mu)
{
  return S1::S1(x, y, c, primes, lpf, mu);
}

int128_t S1(int128_t x,
            int64_t y,
            int64_t c,
            vector<int32_t>& primes,
            FactorTable& factors)
{
  return S1::S1(x, y, c, primes, factors);
}

#endif

} // namespace primecount
