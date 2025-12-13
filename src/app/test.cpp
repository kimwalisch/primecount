///
/// @file  test.cpp
/// @brief primecount integration tests (option: --test).
///        These tests are also used (by the author) for
///        benchmarking code changes.
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <gourdon.hpp>
#include <int128_t.hpp>
#include <PiTable.hpp>
#include <print.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <random>
#include <sstream>

using namespace primecount;

// test: f1(x, threads) == f2(x)
#define TEST0(f1, f2, iters) \
{ \
  std::cout << "Testing " << #f1 << "(x)" << std::flush; \
  int threads = get_num_threads(); \
  int old_percent = -1; \
  int64_t x = 0; \
 \
  /* test small values */ \
  for (x = 0; x < 10000; x++) \
    check_equal(#f1, x, f1 (x, threads), f2 (x)); \
 \
  /* test random increment */ \
  for (int64_t i = 0; i < iters; i++) \
  { \
    check_equal(#f1, x, f1 (x, threads), f2 (x)); \
    int percent = int(100.0 * (i + 1.0) / iters); \
    if (percent > old_percent) \
      std::cout << "\rTesting " << #f1 "(x) " << percent << "%" << std::flush; \
    old_percent = percent; \
    x += dist(gen); \
  } \
 \
  std::cout << std::endl; \
}

// test: f1(x) == f2(x, threads)
#define TEST1(f1, f2, iters) \
{ \
  std::cout << "Testing " << #f1 << "(x)" << std::flush; \
  int threads = get_num_threads(); \
  int old_percent = -1; \
  int64_t x = 0; \
 \
  /* test small values */ \
  for (x = 0; x < 10000; x++) \
    check_equal(#f1, x, f1 (x), f2 (x, threads)); \
 \
  /* test random increment */ \
  for (int64_t i = 0; i < iters; i++) \
  { \
    check_equal(#f1, x, f1 (x), f2 (x, threads)); \
    int percent = int(100.0 * (i + 1.0) / iters); \
    if (percent > old_percent) \
      std::cout << "\rTesting " << #f1 "(x) " << percent << "%" << std::flush; \
    old_percent = percent; \
    x += dist(gen); \
  } \
 \
  std::cout << std::endl; \
}

// test: f1(x, threads) == f2(x, threads)
#define TEST2(f1, f2, iters) \
{ \
  std::cout << "Testing " << #f1 << "(x)" << std::flush; \
  int threads = get_num_threads(); \
  int old_percent = -1; \
  int64_t x = 0; \
 \
  /* test small values */ \
  for (x = 0; x < 10000; x++) \
    check_equal(#f1, x, f1 (x, threads), f2 (x, threads)); \
 \
  /* test random increment */ \
  for (int64_t i = 0; i < iters; i++) \
  { \
    check_equal(#f1, x, f1 (x, threads), f2 (x, threads)); \
    int percent = int(100.0 * (i + 1.0) / iters); \
    if (percent > old_percent) \
      std::cout << "\rTesting " << #f1 "(x) " << percent << "%" << std::flush; \
    old_percent = percent; \
    x += dist(gen); \
  } \
 \
  std::cout << std::endl; \
}

#define TEST_NTH_PRIME(f1, tiny_iters, iters) \
{ \
  std::cout << "Testing " << #f1 "(x) " << std::flush; \
 \
  int64_t n = 1; \
  int64_t prime = 0; \
  int64_t next = tiny_iters; \
  int old_percent = -1; \
 \
  std::random_device rd; \
  std::mt19937 gen(rd()); \
  std::uniform_int_distribution<int64_t> dist(1, 10000000); \
  int threads = get_num_threads(); \
 \
  for (; n < tiny_iters; n++) \
    check_equal(#f1, n, f1 (n, threads), primesieve::nth_prime(n)); \
 \
  for (int64_t i = 0; i < iters; i++) \
  { \
    prime = primesieve::nth_prime(next, prime); \
    check_equal(#f1, n, f1 (n, threads), prime); \
    int percent = int(100.0 * (i + 1.0) / iters); \
    if (percent > old_percent) \
      std::cout << "\rTesting " << #f1 "(x) " << percent << "%" << std::flush; \
    old_percent = percent; \
    next = dist(gen); \
    n += next; \
  } \
 \
  std::cout << std::endl; \
}

namespace {

void check_equal(string_view_t f1,
                 int64_t x,
                 int64_t res1,
                 int64_t res2)
{
  if (res1 != res2)
  {
    std::ostringstream oss;
    oss << f1 << "(" << x << ") = " << res1
        << " is an error, the correct result is " << res2;
    throw primecount_error(oss.str());
  }
}

void test_pi_cache()
{
  std::cout << "Testing pi_cache(x)" << std::flush;

  for (int64_t x = 0; x <= PiTable::max_cached(); x++)
    check_equal("pi_cache", x, pi_cache(x), pi_primesieve(x));

  std::cout << " 100%" << std::endl;
}

} // namespace

namespace primecount {

void test()
{
  set_print(false);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int64_t> dist(1, 10000000);

  try
  {
    test_pi_cache();

    TEST0(pi_legendre,            pi_primesieve,    100);
    TEST2(pi_meissel,             pi_legendre,      500);
    TEST2(pi_lehmer,              pi_meissel,       500);
    TEST1(pi_lmo1,                pi_meissel,        50);
    TEST1(pi_lmo2,                pi_meissel,       200);
    TEST1(pi_lmo3,                pi_meissel,       300);
    TEST1(pi_lmo4,                pi_meissel,       300);
    TEST1(pi_lmo5,                pi_meissel,       600);
    TEST2(pi_lmo_parallel,        pi_meissel,       900);

    TEST2(pi_deleglise_rivat_64,  pi_lmo_parallel, 1500);
#ifdef HAVE_INT128_T
    TEST2(pi_deleglise_rivat_128, pi_lmo_parallel, 1500);
#endif

    TEST2(pi_gourdon_64,          pi_lmo_parallel, 1500);
#ifdef HAVE_INT128_T
    TEST2(pi_gourdon_128,         pi_lmo_parallel, 1500);
#endif

    TEST_NTH_PRIME(nth_prime_64,  10000, 300);
#ifdef HAVE_INT128_T
    TEST_NTH_PRIME(nth_prime_128, 10000, 300);
#endif
  }
  catch (std::exception& e)
  {
    std::cerr << std::endl << e.what() << std::endl;
    std::exit(1);
  }

  std::cout << "All tests passed successfully!" << std::endl;
  std::exit(0);
}

} // namespace
