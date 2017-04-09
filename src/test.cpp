///
/// @file   test.cpp
/// @brief  primecount integration tests.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <primecount.hpp>
#include <primesieve.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <string>
#include <exception>
#include <sstream>
#include <ctime>

#ifdef _OPENMP
  #include <omp.h>
#endif

/// For types: f1(x) , f2(x, threads)
#define CHECK_12(f1, f2) check_equal(#f1, x, f1 (x), f2 (x, get_num_threads()))

/// For types: f1(x, threads) , f2(x, threads)
#define CHECK_22(f1, f2) check_equal(#f1, x, f1 (x, get_num_threads()), f2 (x, get_num_threads()))

#define CHECK_EQUAL(f1, f2, check, iters) \
{ \
  cout << "Testing " << #f1 << "(x)" << flush; \
 \
  /* test for 0 <= x < 10000 */ \
  for (int64_t x = 0; x < 10000; x++) \
    check(f1, f2); \
 \
  int64_t x = 0; \
  /* test using random increment */ \
  for (int64_t i = 0; i < iters; i++, x += get_rand()) \
  { \
    check(f1, f2); \
    double percent = 100.0 * (i + 1.0) / iters; \
    cout << "\rTesting " << #f1 "(x) " << (int) percent << "%" << flush; \
  } \
 \
  cout << endl; \
}

using namespace std;
using namespace primecount;

namespace {

int get_rand()
{
  // 0 <= get_rand() < 10^7
  return (rand() % 10000) * 1000 + 1;
}

void check_equal(const string& f1, int64_t x, int64_t res1, int64_t res2)
{
  if (res1 != res2)
  {
    ostringstream oss;
    oss << f1 << "(" << x << ") = " << res1
        << " is an error, the correct result is " << res2;
    throw primecount_error(oss.str());
  }
}

void check_nth_prime(int64_t iters)
{
  cout << "Testing nth_prime(x)" << flush;

  int64_t n = 0;
  int64_t old_n = 0;
  int64_t old_nth_prime = 0;

  for (; n < 10000; n++)
    check_equal("nth_prime", n, nth_prime(n), primesieve::nth_prime(n));

  // test using random increment
  for (int64_t i = 0; i < iters; i++, n += get_rand())
  {
    int64_t primesieve_nth_prime = primesieve::nth_prime(n - old_n, old_nth_prime);
    check_equal("nth_prime", n, nth_prime(n), primesieve_nth_prime);
    double percent = 100.0 * (i + 1.0) / iters;
    cout << "\rTesting nth_prime(x) " << (int) percent << "%" << flush;

    old_n = n;
    old_nth_prime = primesieve_nth_prime;
  }

  cout << endl;
}

#ifdef _OPENMP

void test_phi_thread_safety(int64_t iters)
{
  cout << "Testing phi(x, a)" << flush;

  int64_t sum1 = 0;
  int64_t sum2 = 0;

  for (int64_t i = 0; i < iters; i++)
    sum1 += pi_legendre(10000000 + i, 1);

  #pragma omp parallel for reduction(+: sum2)
  for (int64_t i = 0; i < iters; i++)
    sum2 += pi_legendre(10000000 + i, 1);

  if (sum1 != sum2)
    throw primecount_error("Error: multi-threaded phi(x, a) is broken.");

  std::cout << "\rTesting phi(x, a) 100%" << endl;
}

#endif

} // namespace

namespace primecount {

bool test()
{
  set_print_status(false); 
  srand(static_cast<unsigned>(time(0)));

  try
  {
#ifdef _OPENMP
    test_phi_thread_safety(100);
#endif

    CHECK_EQUAL(pi_legendre,                  pi_primesieve,      CHECK_22,  100);
    CHECK_EQUAL(pi_meissel,                   pi_legendre,        CHECK_22,  500);
    CHECK_EQUAL(pi_lehmer,                    pi_meissel,         CHECK_22,  500);
    CHECK_EQUAL(pi_lmo1,                      pi_meissel,         CHECK_12,   50);
    CHECK_EQUAL(pi_lmo2,                      pi_meissel,         CHECK_12,  200);
    CHECK_EQUAL(pi_lmo3,                      pi_meissel,         CHECK_12,  300);
    CHECK_EQUAL(pi_lmo4,                      pi_meissel,         CHECK_12,  300);
    CHECK_EQUAL(pi_lmo5,                      pi_meissel,         CHECK_12,  600);
    CHECK_EQUAL(pi_lmo_parallel1,             pi_meissel,         CHECK_22,  400);
    CHECK_EQUAL(pi_lmo_parallel2,             pi_meissel,         CHECK_22,  400);
    CHECK_EQUAL(pi_lmo_parallel3,             pi_meissel,         CHECK_22,  900);
    CHECK_EQUAL(pi_deleglise_rivat1,          pi_lmo_parallel3,   CHECK_12,  600);
    CHECK_EQUAL(pi_deleglise_rivat2,          pi_lmo_parallel3,   CHECK_12,  600);
    CHECK_EQUAL(pi_deleglise_rivat_parallel1, pi_lmo_parallel3,   CHECK_22,  900);
    CHECK_EQUAL(pi_deleglise_rivat_parallel2, pi_lmo_parallel3,   CHECK_22, 1500);

#ifdef HAVE_INT128_T
    CHECK_EQUAL(pi_deleglise_rivat_parallel3, pi_lmo_parallel3,   CHECK_22, 1500);
#endif

    check_nth_prime(300);
  }
  catch (exception& e)
  {
    cerr << endl << e.what() << endl;
    return false;
  }

  cout << "All tests passed successfully!" << endl;
  return true;
}

} // namespace
