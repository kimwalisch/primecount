///
/// @file  primecount-internal.hpp
/// @brief primecount internal functions
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PRIMECOUNT_INTERNAL_HPP
#define PRIMECOUNT_INTERNAL_HPP

#include <int128_t.hpp>
#include <print.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <algorithm>
#include <string>
#include <utility>

namespace primecount {

template<class T>
void unused_param(const T&)
{ }

int64_t pi_lmo1(int64_t x);
int64_t pi_lmo2(int64_t x);
int64_t pi_lmo3(int64_t x);
int64_t pi_lmo4(int64_t x);
int64_t pi_primesieve(int64_t x);

std::string pi(const std::string& x, int threads);
int64_t pi(int64_t x, int threads);
int64_t pi_noprint(int64_t x, int threads);
int64_t pi_deleglise_rivat(int64_t x, int threads);

int64_t pi_cache(int64_t x, bool print = is_print());
int64_t pi_deleglise_rivat_64(int64_t x, int threads, bool print = is_print());
int64_t pi_legendre(int64_t x, int threads, bool print = is_print());
int64_t pi_lehmer(int64_t x, int threads, bool print = is_print());
int64_t pi_lmo5(int64_t x, bool print = is_print());
int64_t pi_lmo_parallel(int64_t x, int threads, bool print = is_print());
int64_t pi_meissel(int64_t x, int threads, bool print = is_print());
int64_t phi(int64_t x, int64_t a, int threads, bool print = is_print());
int64_t P2(int64_t x, int64_t y, int64_t a, int threads, bool print = is_print());
int64_t P3(int64_t x, int64_t y, int64_t a, int threads, bool print = is_print());

int64_t nth_prime(int64_t n, int threads);
int64_t nth_prime_64(int64_t n, int threads);

int64_t Li(int64_t);
int64_t Li_inverse(int64_t);
int64_t RiemannR(int64_t);
int64_t RiemannR_inverse(int64_t);

#ifdef HAVE_INT128_T
  int128_t pi(int128_t x);
  int128_t pi(int128_t x, int threads);
  int128_t pi_deleglise_rivat(int128_t x, int threads);
  int128_t pi_deleglise_rivat_128(int128_t x, int threads, bool print = is_print());
  int128_t P2(int128_t x, int64_t y, int64_t a, int threads, bool print = is_print());

  int128_t nth_prime(int128_t n, int threads);
  int128_t nth_prime_128(int128_t n, int threads);

  int128_t Li(int128_t);
  int128_t Li_inverse(int128_t);
  int128_t RiemannR(int128_t);
  int128_t RiemannR_inverse(int128_t);
#endif

void set_status_precision(int precision);
int get_status_precision(maxint_t x);
void set_alpha(double alpha);
void set_alpha_y(double alpha_y);
void set_alpha_z(double alpha_z);
double get_alpha(maxint_t x, int64_t y);
double get_alpha_y(maxint_t x, int64_t y);
double get_alpha_z(int64_t y, int64_t z);
double get_alpha_lmo(maxint_t x);
double get_alpha_deleglise_rivat(maxint_t x);
std::pair<double, double> get_alpha_gourdon(maxint_t x);
int64_t get_x_star_gourdon(maxint_t x, int64_t y);
maxint_t get_max_x(double alpha_y);
maxint_t to_maxint(const std::string& expr);
double get_time();

} // namespace primecount

namespace {

template <typename T1, typename T2, typename T3>
T2 in_between(T1 min, T2 x, T3 max)
{
  if (x < min || max < min)
    return (T2) min;
  if (x > max)
    return (T2) max;

  return x;
}

inline int ideal_num_threads(int64_t sieve_limit, int threads, int64_t thread_threshold)
{
  thread_threshold = std::max((int64_t) 1, thread_threshold);
  int64_t max_threads = ceil_div(sieve_limit, thread_threshold);
  return in_between(1, threads, max_threads);
}

template <typename T>
double get_percent(T low, T limit)
{
  double percent = (100.0 * low) / std::max((T) 1, limit);
  return in_between(0, percent, 100);
}

template <typename T>
T S2_approx(T x, int64_t pi_y, T p2, T s1)
{
  T pix = primecount::Li(x);
  T s2_approx = pix - s1 - pi_y + 1 + p2;
  s2_approx = std::max(s2_approx, (T) 0);
  return s2_approx;
}

template <typename T>
T D_approx(T x, T sigma, T phi0, T ac, T b)
{
  T pix = primecount::Li(x);
  T d_approx = pix - (ac - b + phi0 + sigma);
  d_approx = std::max(d_approx, (T) 0);
  return d_approx;
}

} // namespace

#endif
