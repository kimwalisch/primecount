///
/// @file  test.cpp
/// @brief primecount integration tests (option: --test).
///        These tests are also used (by the author) for
///        benchmarking code changes.
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <int128_t.hpp>
#include <print.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <random>
#include <sstream>
#include <string>

#ifdef _OPENMP
  #include <omp.h>
#endif

// test: f1(x) == f2(x)
#define TEST0(f1, f2, iters) \
{ \
  cout << "Testing " << #f1 << "(x)" << flush; \
  int64_t x = 0; \
 \
  /* test small values */ \
  for (x = 0; x < 10000; x++) \
    check_equal(#f1, x, f1 (x), f2 (x)); \
 \
  /* test random increment */ \
  for (int64_t i = 0; i < iters; i++) \
  { \
    check_equal(#f1, x, f1 (x), f2 (x)); \
    double percent = 100.0 * (i + 1.0) / iters; \
    cout << "\rTesting " << #f1 "(x) " << (int) percent << "%" << flush; \
    x += dist(gen); \
  } \
 \
  cout << endl; \
}

// test: f1(x) == f2(x, threads)
#define TEST1(f1, f2, iters) \
{ \
  cout << "Testing " << #f1 << "(x)" << flush; \
  int threads = get_num_threads(); \
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
    double percent = 100.0 * (i + 1.0) / iters; \
    cout << "\rTesting " << #f1 "(x) " << (int) percent << "%" << flush; \
    x += dist(gen); \
  } \
 \
  cout << endl; \
}

// test: f1(x, threads) == f2(x, threads)
#define TEST2(f1, f2, iters) \
{ \
  cout << "Testing " << #f1 << "(x)" << flush; \
  int threads = get_num_threads(); \
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
    double percent = 100.0 * (i + 1.0) / iters; \
    cout << "\rTesting " << #f1 "(x) " << (int) percent << "%" << flush; \
    x += dist(gen); \
  } \
 \
  cout << endl; \
}

using namespace std;
using namespace primecount;

namespace {

void check_equal(const string& f1,
                 int64_t x,
                 int64_t res1,
                 int64_t res2)
{
  if (res1 != res2)
  {
    ostringstream oss;
    oss << f1 << "(" << x << ") = " << res1
        << " is an error, the correct result is " << res2;
    throw primecount_error(oss.str());
  }
}

void test_nth_prime(int64_t iters)
{
  cout << "Testing nth_prime(x)" << flush;

  int64_t n = 0;
  int64_t prime = 0;
  int64_t next = 10000;

  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<int64_t> dist(1, 10000000);

  for (; n < next; n++)
    check_equal("nth_prime", n, nth_prime(n), primesieve::nth_prime(n));

  // test random increment
  for (int64_t i = 0; i < iters; i++)
  {
    prime = primesieve::nth_prime(next, prime);
    check_equal("nth_prime", n, nth_prime(n), prime);
    double percent = 100.0 * (i + 1.0) / iters;
    cout << "\rTesting nth_prime(x) " << (int) percent << "%" << flush;
    next = dist(gen);
    n += next;
  }

  cout << endl;
}

#ifdef _OPENMP

void test_phi(int64_t iters)
{
  cout << "Testing phi(x, a)" << flush;

  int64_t sum1 = 0;
  int64_t sum2 = 0;

  #pragma omp parallel for reduction(+: sum1)
  for (int64_t i = 0; i < iters; i++)
    sum1 += pi_legendre(10000000 + i, 1);

  for (int64_t i = 0; i < iters; i++)
    sum2 += pi_legendre(10000000 + i, 1);

  if (sum1 != sum2)
    throw primecount_error("Error: multi-threaded phi(x, a) is broken.");

  cout << "\rTesting phi(x, a) 100%" << endl;
}

#endif

} // namespace

namespace primecount {

void test()
{
  set_print(false);

  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<int64_t> dist(1, 10000000);

  try
  {
#ifdef _OPENMP
    test_phi(100);
#endif

    TEST0(pi_legendre,                  pi_primesieve,    100);
    TEST2(pi_meissel,                   pi_legendre,      500);
    TEST2(pi_lehmer,                    pi_meissel,       500);
    TEST1(pi_lmo1,                      pi_meissel,        50);
    TEST1(pi_lmo2,                      pi_meissel,       200);
    TEST1(pi_lmo3,                      pi_meissel,       300);
    TEST1(pi_lmo4,                      pi_meissel,       300);
    TEST1(pi_lmo5,                      pi_meissel,       600);
    TEST2(pi_lmo_parallel,              pi_meissel,       900);
    TEST1(pi_deleglise_rivat1,          pi_lmo_parallel,  600);
    TEST1(pi_deleglise_rivat2,          pi_lmo_parallel,  600);
    TEST2(pi_deleglise_rivat_parallel1, pi_lmo_parallel, 1500);

#ifdef HAVE_INT128_T
    TEST2(pi_deleglise_rivat_parallel2, pi_lmo_parallel, 1500);
#endif

    test_nth_prime(300);
  }
  catch (exception& e)
  {
    cerr << endl << e.what() << endl;
    exit(1);
  }

  cout << "All tests passed successfully!" << endl;
  exit(0);
}

} // namespace
