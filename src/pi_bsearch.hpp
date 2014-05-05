///
/// @file  pi_bsearch.hpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PI_BSEARCH_HPP
#define PI_BSEARCH_HPP

#include <algorithm>
#include <vector>
#include <cassert>

/// Given a vector with the first n primes and x <= nth prime
/// calculate the number of primes below x using binary search.
/// @pre primes[0] = 0, primes[1] = 2, primes[3] = 3, ...
///
template <typename T1, typename T2>
inline T2 pi_bsearch(const std::vector<T1>& primes, T2 x)
{
  assert(primes[0] == 0);
  // +1 is a correction for primes[0] = 0
  return static_cast<T2>(std::upper_bound(primes.begin() + 1, primes.end(), x) - (primes.begin() + 1));
}

/// Given a vector with the first n primes and x <= nth prime
/// calculate the number of primes below x using binary search.
/// @pre primes[0] = 0, primes[1] = 2, primes[3] = 3, ...
///
template <typename T1, typename T2, typename T3>
inline T3 pi_bsearch(const std::vector<T1>& primes, T2 len, T3 x)
{
  assert(primes[0] == 0);
  // +1 is a correction for primes[0] = 0
  return static_cast<T3>(std::upper_bound(primes.begin() + 1, primes.begin() + len + 1, x) - (primes.begin() + 1));
}

#endif
