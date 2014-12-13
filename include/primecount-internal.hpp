///
/// @file   primecount-internal.hpp
/// @brief  primecount internal function definitions.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PRIMECOUNT_INTERNAL_HPP
#define PRIMECOUNT_INTERNAL_HPP

#include <int128.hpp>
#include <aligned_vector.hpp>
#include <pmath.hpp>

#include <stdint.h>
#include <algorithm>
#include <string>
#include <vector>

namespace primecount {

/// Fastest 64-bit integer type for division.
/// On most Intel CPUs before 2015 unsigned 64-bit division is about
/// 10 percent faster than signed division. It is likely that in a few
/// years years signed and unsigned division will run equally fast.
///
typedef uint64_t intfast64_t;

#if defined(HAVE_INT128_T)

/// Fastest 128-bit integer type for division.
/// On the author's Intel Core-i7 4770 CPU from 2013 using uint128_t
/// instead of int128_t gives 10 percent better performance.
///
typedef uint128_t intfast128_t;

#endif

enum {
  /// Uses all CPU cores.
  MAX_THREADS = -1
};

class PhiCache;

std::string pi(const std::string& x, int threads);

int64_t pi(int64_t x, int threads);

#ifdef HAVE_INT128_T

int128_t pi(int128_t x);

int128_t pi(int128_t x, int threads);

#endif

int64_t pi_deleglise_rivat(int64_t x, int threads);

#ifdef HAVE_INT128_T

int128_t pi_deleglise_rivat(int128_t x);

int128_t pi_deleglise_rivat(int128_t x, int threads);

#endif

int64_t pi_deleglise_rivat1(int64_t x);

int64_t pi_deleglise_rivat2(int64_t x);

#ifdef HAVE_INT128_T

int128_t pi_deleglise_rivat3(int128_t x);

#endif

int64_t pi_deleglise_rivat_parallel1(int64_t x, int threads);

int64_t pi_deleglise_rivat_parallel2(int64_t x, int threads);

#ifdef HAVE_INT128_T

int128_t pi_deleglise_rivat_parallel3(int128_t x, int threads);

#endif

int64_t pi_legendre(int64_t x, int threads);

int64_t pi_lehmer(int64_t x, int threads);

int64_t pi_lehmer2(int64_t x, int threads);

int64_t pi_lmo(int64_t x, int threads);

int64_t pi_lmo1(int64_t x);

int64_t pi_lmo2(int64_t x);

int64_t pi_lmo3(int64_t x);

int64_t pi_lmo4(int64_t x);

int64_t pi_lmo5(int64_t x);

int64_t pi_lmo_parallel1(int64_t x, int threads);

int64_t pi_lmo_parallel2(int64_t x, int threads);

int64_t pi_lmo_parallel3(int64_t x, int threads);

int64_t pi_meissel(int64_t x, int threads);

int64_t pi_primesieve(int64_t x, int threads);

int64_t phi(int64_t x, int64_t a, int threads);

int64_t phi(int64_t x, int64_t a, PhiCache* phiCache);

int64_t Li(int64_t);

int64_t Li_inverse(int64_t);

#ifdef HAVE_INT128_T

int128_t Li(int128_t);

int128_t Li_inverse(int128_t);

#endif

int64_t nth_prime(int64_t n, int threads);

int64_t P2(intfast64_t x, int64_t y, int threads);

#ifdef HAVE_INT128_T

int128_t P2(intfast128_t x, int64_t y, int threads);

#endif

int64_t P2_lehmer(int64_t x, int64_t a, int threads);

int64_t P3(int64_t x, int64_t a, int threads);

double get_wtime();

int validate_threads(int threads);

int validate_threads(int threads, int64_t sieve_limit, int64_t thread_threshold = 100000);

void print_result(const std::string& str, maxint_t res, double time);

void print_seconds(double seconds);

bool print_status();

maxint_t to_maxint(const std::string& expr);

template <typename T>
int get_percent(T low, T limit)
{
  int percent = (int) (100.0 * low / std::max<T>(1, limit));
  return in_between(0, percent, 100);
}

template <typename T>
T S2_approx(T x, int64_t pi_y, T P2, T S1)
{
  T pix = Li(x);
  T S2 = pix - S1 - pi_y + 1 + P2;
  return S2;
}

} // namespace primecount

#endif
