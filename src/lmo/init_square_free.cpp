///
/// @file  init_square_free.cpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <init_square_free.hpp>
#include <pmath.hpp>

#include <stdint.h>
#include <algorithm>
#include <vector>

using namespace std;

namespace primecount {

/// Generate vectors containing n values which satisfy:
/// is_square_free(n) && && !is_prime(n) && primes[i] < least_prime_factor[n].
///
void init_square_free_candidates(vector<vector<int32_t> >* square_free_candidates,
                                 vector<int32_t>& lpf,
                                 vector<int32_t>& mu,
                                 vector<int32_t>& pi,
                                 vector<int32_t>& primes,
                                 int64_t c,
                                 int64_t y)
{
  square_free_candidates->clear();
  square_free_candidates->resize(pi[isqrt(y)], vector<int32_t>(1, 0));

  for (int32_t n = 2; n <= y; n++)
    if (mu[n] != 0 && n != primes[pi[n]])
      for (int32_t i = pi[lpf[n]] - 1; i > c; i--)
        (*square_free_candidates)[i].push_back(n);

#if __cplusplus >= 201103L
  for (size_t i = 0; i < square_free_candidates->size(); i++)
    (*square_free_candidates)[i].shrink_to_fit();
#endif
}

/// Initialize the square free iterators.
/// This version is for use in a single-threaded implementation.
///
void init_square_free_iters(vector<vector<int32_t>::iterator >* iters,
                            vector<vector<int32_t> >& square_free_candidates)
{
  for (int64_t i = 0; i < iters->size(); i++)
    (*iters)[i] = square_free_candidates[i].end() - 1;
}

/// Initialize the square free iterators.
/// This version is for use in a parallel implementation.
///
void init_square_free_iters(vector<vector<int32_t>::iterator >* iters,
                            vector<vector<int32_t> >& square_free_candidates,
                            vector<int32_t>& primes,
                            int64_t c,
                            int64_t x,
                            int64_t y,
                            int64_t low)
{
  for (size_t i = c + 1; i < iters->size(); i++)
  {
    int64_t max_m = x / (primes[i] * low);
    (*iters)[i] = upper_bound(square_free_candidates[i].begin(),
        square_free_candidates[i].end(), max_m) - 1;
  }
}

} // namespace
