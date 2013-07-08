///
/// @file  pi_bsearch.h
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PI_BSEARCH_H
#define PI_BSEARCH_H

#include <algorithm>
#include <iterator>

namespace primecount {

/// Given a std::vector with the first n primes and x <= nth prime
/// calculate the number of primes below x using binary search.
/// Run time: O(log x)
///
template <class ForwardIterator, typename T>
inline T pi_bsearch(ForwardIterator first, ForwardIterator last, const T& x)
{
  return static_cast<T>(
      std::distance(first, std::upper_bound(first, last, x)));
}

} // namespace primecount

#endif
