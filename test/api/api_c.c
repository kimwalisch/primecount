///
/// @file   api_c.c
/// @brief  Test primecount's C API.
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.h>

#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

void check(int OK)
{
  if (OK)
    printf("   OK\n");
  else
  {
    printf("   ERROR\n");
    exit(1);
  }
}

int main(void)
{
  printf("primecount version: %s\n", primecount_version());
  printf("threads: %d\n", primecount_get_num_threads());

  primecount_set_num_threads(3);
  printf("new threads: %d\n", primecount_get_num_threads());

  // Test 64-bit pi(-x)
  int64_t n = -1;
  int64_t res = primecount_pi(n);
  printf("primecount_pi(%"PRId64") = %"PRId64, n, res);
  check(res == 0);

  n = -9223372036854775807;
  res = primecount_pi(n);
  printf("primecount_pi(%"PRId64") = %"PRId64, n, res);
  check(res == 0);

  n = (int64_t) 1e10;
  res = primecount_pi(n);
  printf("primecount_pi(%"PRId64") = %"PRId64, n, res);
  check(res == 455052511);

  pc_int128_t n128;
  n128.lo = (uint64_t) 1e9;
  n128.hi = 0;
  pc_int128_t res128 = primecount_pi_128(n128);
  printf("primecount_pi_128(1e9) = %"PRId64, res128.lo);
  check(res128.lo == 50847534 && res128.hi == 0);

  // Check x >= primecount max x of 10^31.
  // primecount must detect issue and return -1 error.
  n128.lo = 0;
  n128.hi = 1ull << 50;
  res128 = primecount_pi_128(n128);
  printf("primecount_pi_128(2^114) = %"PRId64, res128.lo);
  check(res128.hi == -1 && ~res128.lo == 0);

  n = 455052511;
  res = primecount_nth_prime(n);
  printf("primecount_nth_prime(%"PRId64") = %"PRId64, n, res);
  check(res == 9999999967);

  // nth_prime(-1) is an error and should hence return -1
  // which indicates an error in the libprimecount C API.
  n = -1;
  res = primecount_nth_prime(n);
  printf("primecount_nth_prime(%"PRId64") = %"PRId64, n, res);
  check(res == -1);

  n128.lo = (uint64_t) 1e9;
  n128.hi = 0;
  res128 = primecount_nth_prime_128(n128);
  printf("primecount_nth_prime_128(1e9) = %"PRId64, res128.lo);
  check(res128.lo == 22801763489 && res128.hi == 0);

  // Check n >= primecount max n of ~ 10^29.
  // primecount must detect issue and return -1 error.
  n128.lo = 0;
  n128.hi = 1ull << 50;
  res128 = primecount_nth_prime_128(n128);
  printf("primecount_nth_prime_128(1e9) = %"PRId64, res128.lo);
  check(res128.hi == -1 && ~res128.lo == 0);

  n = (int64_t) 1e12;
  int64_t a = 78498;
  res = primecount_phi(n, a);
  printf("primecount_phi(%"PRId64", %"PRId64") = %"PRId64, n , a, res);
  check(res == 37607833521);

  n = -1;
  res = primecount_phi(n, a);
  printf("primecount_phi(%"PRId64", %"PRId64") = %"PRId64, n , a, res);
  check(res == 0);

  printf("\n");
  printf("All tests passed successfully!\n");

  return 0;
}
