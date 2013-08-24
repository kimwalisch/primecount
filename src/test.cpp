///
/// @file  test.cpp
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.h>
#include <primesieve/soe/ParallelPrimeSieve.h>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <string>
#include <ctime>

using namespace std;

namespace primecount {

void assert_equal(const std::string& f1_name, int64_t x, int64_t f1_res, int64_t f2_res)
{
  if (f1_res != f2_res)
  {
    cerr << endl << f1_name << "(" << x << ") = "  << f1_res
         << " is an error, the correct result is " << f2_res << std::endl;
    exit(1);
  }
}

/// 0 <= get_rand() < 10^7
int get_rand()
{
  return (rand() % 10000) * 1000 + 1;
}

template <typename F>
void check_for_equality(const std::string& f1_name, F f1, F f2, int64_t iters)
{
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
  ParallelPrimeSieve pps;
  int64_t prime = pps.nthPrime(x);
  return prime;
}

void test()
{
  srand(static_cast<unsigned int>(time(0)));

  check_for_equality("pi_legendre",   pi_legendre,   pi_primesieve, 100);
  check_for_equality("pi_meissel",    pi_meissel,    pi_legendre,   500);
  check_for_equality("pi_lehmer",     pi_lehmer,     pi_meissel,    500);
  check_for_equality("pi_lmo_simple", pi_lmo_simple, pi_meissel,    200);
  check_for_equality("nth_prime",     nth_prime,     pps_nth_prime, 100);

  cout << "All tests passed successfully!" << endl;
  exit(0);
}

} // namespace primecount
