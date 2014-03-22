///
/// @file  test.cpp
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
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

using namespace std;
using primecount::MAX_THREADS;

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

template <typename F>
void check_equal(const string& f1_name, F f1, F f2, int64_t iters)
{
  srand(static_cast<unsigned int>(time(0)));
  cout << "Testing " << (f1_name + "(x)") << flush;

  // test for 0 <= x < iters
  for (int64_t x = 0; x < iters; x++)
    assert_equal(f1_name, x, f1(x, MAX_THREADS), f2(x, MAX_THREADS));

  int64_t x = 0;
  // test using random increment
  for (int64_t i = 0; i < iters; i++, x += get_rand())
    assert_equal(f1_name, x, f1(x, MAX_THREADS), f2(x, MAX_THREADS));

  cout << " correct" << endl;
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
  try
  {
    check_equal("pi_legendre", pi_legendre, pi_primesieve, 100);
    check_equal("pi_meissel",  pi_meissel,  pi_legendre,   500);
    check_equal("pi_lehmer",   pi_lehmer,   pi_meissel,    500);
    check_equal("pi_lmo1",     pi_lmo1,     pi_lehmer,     500);
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
