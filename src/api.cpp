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
#include <sstream>
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

// Below 10^7 LMO is faster than Gourdon's algorithm
const int lmo_threshold_ = 10000000;

} // namespace

namespace primecount {

int64_t pi(int64_t x)
{
  return pi(x, get_num_threads());
}

int64_t pi(int64_t x, int threads)
{
  if (x <= lmo_threshold_)
    return pi_lmo5(x);
  else
  {
#ifdef ENABLE_MPI
    // So far only the Deleglise-Rivat algorithm has been distributed
    if (mpi_num_procs() > 1)
      return pi_deleglise_rivat_64(x, threads);
#endif

    return pi_gourdon_64(x, threads);
  }
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
  else
  {
#ifdef ENABLE_MPI
    // So far only the Deleglise-Rivat algorithm has been distributed
    if (mpi_num_procs() > 1)
      return pi_deleglise_rivat_128(x, threads);
#endif

    return pi_gourdon_128(x, threads);
  }
}

#endif

string pi(const string& x)
{
  return pi(x, get_num_threads());
}

string pi(const string& x, int threads)
{
  maxint_t pi_x = pi(to_maxint(x), threads);
  ostringstream oss;
  oss << pi_x;
  return oss.str();
}

int64_t pi_deleglise_rivat(int64_t x, int threads)
{
  return pi_deleglise_rivat_64(x, threads);
}

int64_t pi_gourdon(int64_t x, int threads)
{
  return pi_gourdon_64(x, threads);
}

#ifdef HAVE_INT128_T

int128_t pi_deleglise_rivat(int128_t x, int threads)
{
  // use 64-bit if possible
  if (x <= numeric_limits<int64_t>::max())
    return pi_deleglise_rivat_64((int64_t) x, threads);
  else
    return pi_deleglise_rivat_128(x, threads);
}

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
  double max_x = pow(pow(2.0, 62.0) * alpha_y, 3.0 / 2.0);
  return (int128_t) max_x; 
#else
  unused_param(alpha_y); 
  return numeric_limits<int64_t>::max();
#endif
}

std::string get_max_x()
{
  ostringstream oss;
  oss << get_max_x(1.0);
  return oss.str();
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
