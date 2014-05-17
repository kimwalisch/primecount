///
/// @file  test.cpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "internal.hpp"

#include <primecount.hpp>
#include <primesieve.hpp>

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

using namespace std;

namespace {

int get_rand()
{
  // 0 <= get_rand() < 10^7
  return (rand() % 10000) * 1000 + 1;
}

void assert_equal(const string& f1, int64_t x, int64_t res1, int64_t res2)
{
  if (res1 != res2)
  {
    ostringstream oss;
    oss << f1 << "(" << x << ") = " << res1
        << " is an error, the correct result is " << res2;
    throw runtime_error(oss.str());
  }
}

#define ASSERT_EQUAL(f1, f2, iters) \
{ \
  cout << "Testing " << #f1 "(x)" << flush; \
 \
  /* test for 0 <= x < 1000 */ \
  for (int64_t x = 0; x < 1000; x++) \
    assert_equal(#f1, x, f1 (x), f2 (x)); \
 \
  int64_t x = 0; \
  /* test using random increment */ \
  for (int64_t i = 0; i < iters; i++, x += get_rand()) \
  { \
    assert_equal(#f1, x, f1 (x), f2 (x)); \
    double percent = 100.0 * (i + 1.0) / iters; \
    cout << "\rTesting " << #f1 "(x) " << static_cast<int>(percent) << "%" << flush; \
  } \
 \
  cout << endl; \
}

#ifdef _OPENMP

void test_phi_thread_safety(int64_t iters)
{
  cout << "Testing phi(x, a)" << flush;

  int nested_threads = 2;
  int64_t single_thread_sum = 0;
  int64_t multi_thread_sum = 0;
  int64_t base = 1000000;

  omp_set_nested(true);

  #pragma omp parallel for reduction(+: multi_thread_sum)
  for (int64_t i = 0; i < iters; i++)
    multi_thread_sum += primecount::pi_legendre(base + i, nested_threads);

  omp_set_nested(false);

  for (int64_t i = 0; i < iters; i++)
    single_thread_sum += primecount::pi_legendre(base + i);

  if (multi_thread_sum != single_thread_sum)
    throw runtime_error("Error: multi-threaded phi(x, a) is broken.");

  std::cout << "\rTesting phi(x, a) 100%" << endl;
}

#endif

} // namespace

namespace primecount {

bool test()
{
  srand(static_cast<unsigned int>(time(0)));
  try
  {
#ifdef _OPENMP
    test_phi_thread_safety(100);
#endif
    ASSERT_EQUAL(pi_legendre,      pi_primesieve, 100);
    ASSERT_EQUAL(pi_meissel,       pi_legendre,   400);
    ASSERT_EQUAL(pi_lehmer,        pi_meissel,    400);
    ASSERT_EQUAL(pi_lehmer2,       pi_lehmer,     200);
    ASSERT_EQUAL(pi_lmo1,          pi_meissel,    400);
    ASSERT_EQUAL(pi_lmo2,          pi_meissel,    200);
    ASSERT_EQUAL(pi_lmo3,          pi_meissel,    400);
    ASSERT_EQUAL(pi_lmo4,          pi_meissel,    400);
    ASSERT_EQUAL(pi_lmo5,          pi_meissel,    400);
    ASSERT_EQUAL(pi_lmo_parallel1, pi_meissel,    400);
    ASSERT_EQUAL(pi_lmo_parallel2, pi_meissel,    400);
    ASSERT_EQUAL(nth_prime,        primesieve::parallel_nth_prime, 100);
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
