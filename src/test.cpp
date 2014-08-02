///
/// @file   test.cpp
/// @brief  primecount integration tests.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <ptypes.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <sstream>
#include <ctime>

#ifdef _OPENMP
  #include <omp.h>
#endif

/// For types: f1(x) , f2(x)
#define CHECK_11(f1, f2) check_equal(#f1, x, f1 (x), f2 (x))

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
using primesieve::parallel_nth_prime;

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
    throw runtime_error(oss.str());
  }
}

void test_phi_thread_safety(int64_t iters)
{
#ifdef _OPENMP
  cout << "Testing phi(x, a)" << flush;

  int nested_threads = 2;
  int64_t single_thread_sum = 0;
  int64_t multi_thread_sum = 0;
  int64_t base = 1000000;

  omp_set_nested(true);

  #pragma omp parallel for reduction(+: multi_thread_sum)
  for (int64_t i = 0; i < iters; i++)
    multi_thread_sum += pi_legendre(base + i, nested_threads);

  omp_set_nested(false);

  for (int64_t i = 0; i < iters; i++)
    single_thread_sum += pi_legendre(base + i, 1);

  if (multi_thread_sum != single_thread_sum)
    throw runtime_error("Error: multi-threaded phi(x, a) is broken.");

  std::cout << "\rTesting phi(x, a) 100%" << endl;
#endif
}

} // namespace

namespace primecount {

bool test()
{
  srand((unsigned) time(0));
  try
  {
    test_phi_thread_safety(100);

    CHECK_EQUAL(pi_legendre,                  pi_primesieve,      CHECK_22, 100);
    CHECK_EQUAL(pi_meissel,                   pi_legendre,        CHECK_22, 400);
    CHECK_EQUAL(pi_lehmer,                    pi_meissel,         CHECK_22, 400);
    CHECK_EQUAL(pi_lehmer2,                   pi_lehmer,          CHECK_22, 200);
    CHECK_EQUAL(pi_lmo1,                      pi_meissel,         CHECK_12, 200);
    CHECK_EQUAL(pi_lmo2,                      pi_meissel,         CHECK_12, 200);
    CHECK_EQUAL(pi_lmo3,                      pi_meissel,         CHECK_12, 300);
    CHECK_EQUAL(pi_lmo4,                      pi_meissel,         CHECK_12, 300);
    CHECK_EQUAL(pi_lmo5,                      pi_meissel,         CHECK_12, 400);
    CHECK_EQUAL(pi_lmo_parallel1,             pi_meissel,         CHECK_22, 400);
    CHECK_EQUAL(pi_lmo_parallel2,             pi_meissel,         CHECK_22, 400);
    CHECK_EQUAL(pi_lmo_parallel3,             pi_meissel,         CHECK_22, 400);
    CHECK_EQUAL(pi_deleglise_rivat1,          pi_lmo_parallel3,   CHECK_12, 600);
    CHECK_EQUAL(pi_deleglise_rivat2,          pi_lmo_parallel3,   CHECK_12, 600);
    CHECK_EQUAL(pi_deleglise_rivat3,          pi_lmo_parallel3,   CHECK_12, 600);
#ifdef HAVE_INT128_T
    CHECK_EQUAL(pi_deleglise_rivat4,          pi_lmo_parallel3,   CHECK_12, 600);
#endif
    CHECK_EQUAL(pi_deleglise_rivat_parallel1, pi_lmo_parallel3,   CHECK_22, 900);
    CHECK_EQUAL(pi_deleglise_rivat_parallel2, pi_lmo_parallel3,   CHECK_22, 900);
    CHECK_EQUAL(pi_deleglise_rivat_parallel3, pi_lmo_parallel3,   CHECK_22, 900);
#ifdef HAVE_INT128_T
    CHECK_EQUAL(pi_deleglise_rivat_parallel4, pi_lmo_parallel3,   CHECK_22, 900);
#endif
    CHECK_EQUAL(nth_prime,                    parallel_nth_prime, CHECK_11,  70);
  }
  catch (runtime_error& e)
  {
    cerr << endl << e.what() << endl;
    return false;
  }

  cout << "All tests passed successfully!" << endl;
  return true;
}

} // namespace primecount
