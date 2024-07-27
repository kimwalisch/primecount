///
/// @file  P3.cpp
/// @brief P3(x, a) is the 3rd partial sieve function, it is used
///        in Lehmer's prime counting formula.
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <generate_primes.hpp>
#include <imath.hpp>
#include <macros.hpp>
#include <PiTable.hpp>
#include <print.hpp>

#include <stdint.h>

namespace primecount {

/// P3(x, a) counts the numbers <= x that have exactly 3
/// prime factors each exceeding the a-th prime.
/// Memory usage: O(sqrt(x))
///
int64_t P3(int64_t x,
           int64_t y,
           int64_t a,
           int threads,
           bool is_print)
{
  ASSERT(a == pi_noprint(y, threads));
  double time;

  if (is_print)
  {
    print("");
    print("=== P3(x, a) ===");
    time = get_time();
  }

  int64_t sum = 0;
  int64_t x13 = iroot<3>(x);

  if (y <= x13)
  {
    int64_t max_prime = std::max(x13, isqrt(x / y));
    int64_t max_pix = std::max(x13, x / (y * y));
    auto primes = generate_primes<int32_t>(max_prime);
    PiTable pi(max_pix, threads);
    int64_t pi_x13 = pi[x13];

    // These load balancing settings work well on my
    // dual-socket AMD EPYC 7642 server with 192 CPU cores.
    int64_t thread_threshold = 100;
    threads = ideal_num_threads(pi_x13, threads, thread_threshold);

    #pragma omp parallel for schedule(dynamic, 16) num_threads(threads) reduction(+: sum)
    for (int64_t i = a + 1; i <= pi_x13; i++)
    {
      int64_t xi = x / primes[i];
      int64_t bi = pi[isqrt(xi)];

      for (int64_t j = i; j <= bi; j++)
        sum += pi[xi / primes[j]] - (j - 1);
    }
  }

  if (is_print)
    print("P3", sum, time);

  return sum;
}

} // namespace
