///
/// @file  test.cpp
/// @brief primecount integration tests (option: --test).
///        These tests are also used (by the author) for
///        benchmarking code changes.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <gourdon.hpp>
#include <int128_t.hpp>
#include <macros.hpp>
#include <PiTable.hpp>
#include <print.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <random>

using namespace primecount;

namespace {

NOINLINE void check_equal(string_view_t f1,
                          int64_t x,
                          int64_t res1,
                          int64_t res2)
{
  if (res1 != res2)
  {
    std::cerr << std::endl << f1 << "(" << x << ") = " << res1
              << " is an error, the correct result = " << res2
              << std::endl;
    std::exit(1);
  }
}

void test_pi_cache()
{
  std::cout << "Testing pi_cache(x)" << std::flush;

  for (int64_t x = 0; x <= PiTable::max_cached(); x++)
    check_equal("pi_cache", x, pi_cache(x), pi_primesieve(x));

  std::cout << " 100%" << std::endl;
}

/// test: f1(x) == f2(x)
template <typename F1, typename F2>
void test_pi(string_view_t name,
             F1 f1,
             F2 f2,
             int64_t iters,
             std::mt19937& gen,
             std::uniform_int_distribution<int64_t>& dist)
{
  std::cout << "Testing " << name << "(x)" << std::flush;
  int old_percent = -1;
  int64_t x = 0;

  // test small values
  for (x = 0; x < 10000; x++)
    check_equal(name, x, f1(x), f2(x));

  // test random increment
  for (int64_t i = 0; i < iters; i++)
  {
    check_equal(name, x, f1(x), f2(x));
    int percent = int(100.0 * (i + 1.0) / iters);
    if (percent > old_percent)
      std::cout << "\rTesting " << name << "(x) " << percent << "%" << std::flush;
    old_percent = percent;
    x += dist(gen);
  }

  std::cout << std::endl;
}

/// test: f(n) == primesieve::nth_prime(n)
template <typename F>
void test_nth_prime(string_view_t name,
                    F f,
                    int64_t tiny_iters,
                    int64_t iters,
                    std::mt19937& gen,
                    std::uniform_int_distribution<int64_t>& dist)
{
  std::cout << "Testing " << name << "(x) " << std::flush;
  int old_percent = -1;
  int64_t n = 1;
  int64_t prime = 0;
  int64_t next = tiny_iters;

  for (; n < tiny_iters; n++)
    check_equal(name, n, f(n), primesieve::nth_prime(n));

  for (int64_t i = 0; i < iters; i++)
  {
    prime = primesieve::nth_prime(next, prime);
    check_equal(name, n, f(n), prime);
    int percent = int(100.0 * (i + 1.0) / iters);
    if (percent > old_percent)
      std::cout << "\rTesting " << name << "(x) " << percent << "%" << std::flush;
    old_percent = percent;
    next = dist(gen);
    n += next;
  }

  std::cout << std::endl;
}

} // namespace

// Each call site specifies the arguments of the tested functions,
// e.g. TEST(pi_lmo2, (x, alpha), pi_meissel, (x, threads), 200)
#define TEST(f1, args1, f2, args2, iters) \
  test_pi(#f1, [&](int64_t x) { return f1 args1; }, \
               [&](int64_t x) { return f2 args2; }, \
          iters, gen, dist)

#define TEST_NTH_PRIME(f1, args1, tiny_iters, iters) \
  test_nth_prime(#f1, [&](int64_t n) { return f1 args1; }, \
                 tiny_iters, iters, gen, dist)

namespace primecount {

/// Run all tests and exit.
/// Calls std::exit(1) if any error occurs.
///
void test(double alpha)
{
  set_print(false);
  int threads = get_num_threads();

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int64_t> dist(1, 10000000);

  try
  {
    test_pi_cache();

    TEST(pi_legendre,            (x, threads),        pi_primesieve,   (x),           100);
    TEST(pi_meissel,             (x, threads),        pi_legendre,     (x, threads),  500);
    TEST(pi_lehmer,              (x, threads),        pi_meissel,      (x, threads),  500);
    TEST(pi_lmo1,                (x),                 pi_meissel,      (x, threads),   50);
    TEST(pi_lmo2,                (x, alpha),          pi_meissel,      (x, threads),  200);
    TEST(pi_lmo3,                (x, alpha),          pi_meissel,      (x, threads),  300);
    TEST(pi_lmo4,                (x, alpha),          pi_meissel,      (x, threads),  300);
    TEST(pi_lmo5,                (x, alpha),          pi_meissel,      (x, threads),  600);
    TEST(pi_lmo_parallel,        (x, threads, alpha), pi_meissel,      (x, threads),  900);

    TEST(pi_deleglise_rivat_64,  (x, threads, alpha), pi_lmo_parallel, (x, threads), 1500);
#ifdef HAVE_INT128_T
    TEST(pi_deleglise_rivat_128, (x, threads, alpha), pi_lmo_parallel, (x, threads), 1500);
#endif

    TEST(pi_gourdon_64,          (x, threads, alpha), pi_lmo_parallel, (x, threads), 1500);
#ifdef HAVE_INT128_T
    TEST(pi_gourdon_128,         (x, threads, alpha), pi_lmo_parallel, (x, threads), 1500);
#endif

    TEST_NTH_PRIME(nth_prime_64,  (n, threads, alpha), 10000, 300);
#ifdef HAVE_INT128_T
    TEST_NTH_PRIME(nth_prime_128, (n, threads, alpha), 10000, 300);
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
