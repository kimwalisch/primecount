///
/// @file   alpha_z.cpp
/// @brief  Test the alpha_z tuning factor in Gourdon's algorithm.
///         z = y * alpha_z
///         By computing pi(x) using different alpha_z tuning
///         factors we can make sure that all array sizes
///         (and other bounds) are accurate.
///
/// Copyright (C) 2023 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <gourdon.hpp>
#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>

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

  // Test small values of x
  {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int64_t> dist(100, 1000);

    for (int i = 0; i < 100; i++)
    {
      int64_t x = dist(gen);
      int64_t res1 = pi_cache(x);

      for (double alpha_z = 1; alpha_z <= iroot<6>(x); alpha_z++)
      {
        set_alpha_z(alpha_z);
        int64_t res2 = pi_gourdon_64(x, threads);

        std::cout << "pi_gourdon_64(" << x << ") = " << res2;
        check(res1 == res2);
      }
    }

  #ifdef HAVE_INT128_T
    for (int i = 0; i < 100; i++)
    {
      int64_t x = dist(gen);
      int64_t res1 = pi_cache(x);

      for (double alpha_z = 1; alpha_z <= iroot<6>(x); alpha_z++)
      {
        set_alpha_z(alpha_z);
        int64_t res2 = pi_gourdon_128(x, threads);

        std::cout << "pi_gourdon_128(" << x << ") = " << res2;
        check(res1 == res2);
      }
    }
  #endif
  }

  // Test medium values of x
  {
    int64_t min = (int64_t) 1e8;
    int64_t max = min * 2;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int64_t> dist(min, max);

    for (int i = 0; i < 50; i++)
    {
      int64_t x = dist(gen);
      int64_t res1 = pi_meissel(x, threads);

      for (double alpha_z = 1; alpha_z <= iroot<6>(x); alpha_z++)
      {
        set_alpha_z(alpha_z);
        int64_t res2 = pi_gourdon_64(x, threads);

        std::cout << "pi_gourdon_64(" << x << ") = " << res2;
        check(res1 == res2);
      }
    }

  #ifdef HAVE_INT128_T
    for (int i = 0; i < 50; i++)
    {
      int64_t x = dist(gen);
      int64_t res1 = pi_meissel(x, threads);

      for (double alpha_z = 1; alpha_z <= iroot<6>(x); alpha_z++)
      {
        set_alpha_z(alpha_z);
        int64_t res2 = pi_gourdon_128(x, threads);

        std::cout << "pi_gourdon_128(" << x << ") = " << res2;
        check(res1 == res2);
      }
    }
  #endif
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
