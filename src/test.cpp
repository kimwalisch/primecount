///
/// @file  test.cpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

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
using namespace primecount;

namespace {

void assert_equal(const string& f1_name, int64_t x, int64_t f1_res, int64_t f2_res)
{
  if (f1_res != f2_res)
  {
    ostringstream oss;
    oss << f1_name << "(" << x << ") = " << f1_res
        << " is an error, the correct result is " << f2_res;
    throw runtime_error(oss.str());
  }
}

/// 0 <= get_rand() < 10^7
int get_rand()
{
  return (rand() % 10000) * 1000 + 1;
}

#ifdef _OPENMP

void check_phi_thread_safety(int64_t iters)
{
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
    single_thread_sum += pi_legendre(base + i);

  if (multi_thread_sum != single_thread_sum)
    throw runtime_error("Error: multi-threaded phi(x, a) is broken.");

  std::cout << "\rTesting phi(x, a) 100%" << endl;
}

#endif

template <typename F>
void check_equal(const string& f1_name, F f1, F f2, int64_t iters)
{
  cout << "Testing " << (f1_name + "(x)") << flush;

  // test for 0 <= x < 1000
  for (int64_t x = 0; x < 1000; x++)
    assert_equal(f1_name, x, f1(x, MAX_THREADS), f2(x, MAX_THREADS));

  int64_t x = 0;
  // test using random increment
  for (int64_t i = 0; i < iters; i++, x += get_rand())
  {
    assert_equal(f1_name, x, f1(x, MAX_THREADS), f2(x, MAX_THREADS));
    double percent = 100.0 * (i + 1.0) / iters;
    cout << "\rTesting " << (f1_name + "(x) ") << static_cast<int>(percent) << "%" << flush;
  }

  cout << endl;
}

int64_t pps_nth_prime(int64_t x, int)
{
  int64_t prime = primesieve::parallel_nth_prime(x);
  return prime;
}

} // namespace

namespace primecount {

bool test()
{
  srand(static_cast<unsigned int>(time(0)));
  try
  {
#ifdef _OPENMP
    check_phi_thread_safety(100);
#endif
    check_equal("pi_legendre", pi_legendre, pi_primesieve, 100);
    check_equal("pi_meissel",  pi_meissel,  pi_legendre,   400);
    check_equal("pi_lehmer",   pi_lehmer,   pi_meissel,    400);
    check_equal("pi_lehmer2",  pi_lehmer2,  pi_lehmer,     200);
    check_equal("pi_lmo1",     pi_lmo1,     pi_lehmer,     400);
    check_equal("pi_lmo2",     pi_lmo2,     pi_lehmer,     200);
    check_equal("pi_lmo3",     pi_lmo3,     pi_lehmer,     400);
    check_equal("pi_lmo4",     pi_lmo4,     pi_lehmer,     400);
    check_equal("pi_lmo5",     pi_lmo5,     pi_lehmer,     400);
    check_equal("nth_prime",   nth_prime,   pps_nth_prime, 100);
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
