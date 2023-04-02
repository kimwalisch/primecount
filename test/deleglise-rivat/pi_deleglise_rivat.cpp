///
/// @file   pi_deleglise_rivat.cpp
/// @brief  Test the pi_deleglise_rivat_64(x) and
///         pi_deleglise_rivat_128(x) functions.
///
/// Copyright (C) 2023 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <PiTable.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <random>

using namespace primecount;

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  int threads = get_num_threads();

  // Test negative x
  {
    int64_t x = -1;
    int64_t res = pi_deleglise_rivat_64(x, threads);
    std::cout << "pi_deleglise_rivat_64(" << x << ") = " << res;
    check(res == 0);
  }

#ifdef HAVE_INT128_T
  {
    int128_t x = -1;
    int128_t res = pi_deleglise_rivat_128(x, threads);
    std::cout << "pi_deleglise_rivat_128(" << x << ") = " << res;
    check(res == 0);

    // Test if cast in pi_deleglise_rivat(x) supports x <= -2^64
    x = -1 * (((int128_t) 1) << 100);
    res = pi_deleglise_rivat(x, threads);
    std::cout << "pi_deleglise_rivat(" << x << ") = " << res;
    check(res == 0);
  }
#endif

  // Test small x
  for (int64_t x = 0; x <= PiTable::max_cached(); x++)
  {
    int64_t res1 = pi_cache(x);
    int64_t res2 = pi_deleglise_rivat_64(x, threads);
    std::cout << "pi_deleglise_rivat_64(" << x << ") = " << res2;
    check(res2 == res1);

    #ifdef HAVE_INT128_T
      int128_t res3 = pi_deleglise_rivat_128(x, threads);
      std::cout << "pi_deleglise_rivat_128(" << x << ") = " << res3;
      check(res3 == res1);
    #endif
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int64_t> dist(0, 1 << 27);

  // Test medium x
  for (int i = 0; i < 1000; i++)
  {
    int64_t x = dist(gen);
    int64_t res1 = pi_meissel(x, threads);
    int64_t res2 = pi_deleglise_rivat_64(x, threads);
    std::cout << "pi_deleglise_rivat_64(" << x << ") = " << res2;
    check(res2 == res1);

    #ifdef HAVE_INT128_T
      int128_t res3 = pi_deleglise_rivat_128(x, threads);
      std::cout << "pi_deleglise_rivat_128(" << x << ") = " << res3;
      check(res3 == res1);
    #endif
  }

  {
    // Test larger computation: pi(1e11)
    int64_t x = 100000000000ll;
    int64_t res = pi_deleglise_rivat_64(x, threads);
    std::cout << "pi_deleglise_rivat_64(" << x << ") = " << res;
    check(res == 4118054813ll);
  }

#ifdef HAVE_INT128_T
  {
    // Test larger computation: pi(1e12)
    int128_t x = 1000000000000ll;
    int128_t res = pi_deleglise_rivat_128(x, threads);
    std::cout << "pi_deleglise_rivat_128(" << x << ") = " << res;
    check(res == 37607912018ll);
  }
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
