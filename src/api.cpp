///
/// @file  api.cpp
///        primecount's C++ API.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <gourdon.hpp>
#include <int128_t.hpp>

#include <cmath>
#include <limits>
#include <string>
#include <stdint.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

namespace {

#ifdef _OPENMP
  int threads_ = 0;
#endif

} // namespace

namespace primecount {

int64_t pi(int64_t x)
{
  return pi(x, get_num_threads());
}

int64_t pi(int64_t x, int threads)
{
  return pi_gourdon_64(x, threads);
}

#ifdef HAVE_INT128_T

int128_t pi(int128_t x)
{
  return pi(x, get_num_threads());
}

int128_t pi(int128_t x, int threads)
{
  // use 64-bit if possible
  if (x <= numeric_limits<int64_t>::max())
    return pi((int64_t) x, threads);

  return pi_gourdon_128(x, threads);
}

#endif

string pi(const string& x)
{
  return pi(x, get_num_threads());
}

string pi(const string& x, int threads)
{
  maxint_t n = to_maxint(x);
  maxint_t res = pi(n, threads);
  return to_str(res);
}

int64_t pi_gourdon(int64_t x, int threads)
{
  return pi_gourdon_64(x, threads);
}

#ifdef HAVE_INT128_T

int128_t pi_gourdon(int128_t x, int threads)
{
  // use 64-bit if possible
  if (x <= numeric_limits<int64_t>::max())
    return pi_gourdon_64((int64_t) x, threads);
  else
    return pi_gourdon_128(x, threads);
}

#endif

int64_t nth_prime(int64_t n)
{
  return nth_prime(n, get_num_threads());
}

int64_t phi(int64_t x, int64_t a)
{
  return phi(x, a, get_num_threads());
}

string primecount_version()
{
  return PRIMECOUNT_VERSION;
}

/// Returns the largest x supported by pi(x).
/// The S2_hard, P2, B and D functions are limited by:
/// x / y <= 2^62, with y = x^(1/3) * alpha_y
/// Hence x^(2/3) / alpha_y <= 2^62
/// x <= (2^62 * alpha_y)^(3/2)
///
maxint_t get_max_x(double alpha_y)
{
#ifdef HAVE_INT128_T
  double max_x = pow((1ull << 62) * alpha_y, 3.0 / 2.0);
  return (int128_t) max_x; 
#else
  unused_param(alpha_y); 
  return numeric_limits<int64_t>::max();
#endif
}

std::string get_max_x()
{
#ifdef HAVE_INT128_T
  // 10^31
  return "10000000000000000000000000000000";
#else
  // 2^63-1
  return "9223372036854775807";
#endif
}

int get_num_threads()
{
#ifdef _OPENMP
  if (threads_)
    return threads_;
  else
    return max(1, omp_get_max_threads());
#else
  return 1;
#endif
}

void set_num_threads(int threads)
{
#ifdef _OPENMP
  threads_ = in_between(1, threads, omp_get_max_threads());
#endif
  primesieve::set_num_threads(threads);
}

} // namespace
