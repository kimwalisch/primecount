/*
 * @file  primecount.hpp
 * @brief primecount C API. primecount is a C/C++ library for
 *        counting the number of primes <= x using highly
 *        optimized implementations of the combinatorial type
 *        prime counting function algorithms.
 *
 * Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
 *
 * This file is distributed under the BSD License.
 */

#ifndef PRIMECOUNT_H
#define PRIMECOUNT_H

#include <stdint.h>

#define PRIMECOUNT_VERSION "5.4"
#define PRIMECOUNT_VERSION_MAJOR 5
#define PRIMECOUNT_VERSION_MINOR 4

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Count the number of primes <= x using Xavier Gourdon's
 * algorithm. Uses all CPU cores by default.
 * Returns -1 if an error occurs.
 * 
 * Run time: O(x^(2/3) / (log x)^2)
 * Memory usage: O(x^(1/3) * (log x)^3)
 */
int64_t primecount_pi(int64_t x);

/*
 * Partial sieve function (a.k.a. Legendre-sum).
 * phi(x, a) counts the numbers <= x that are not divisible
 * by any of the first a primes.
 * Returns -1 if an error occurs.
 */
int64_t primecount_phi(int64_t x, int64_t a);

/*
 * Find the nth prime using a combination of the prime counting
 * function and the sieve of Eratosthenes.
 * @pre n <= 216289611853439384
 * Returns -1 if an error occurs.
 * 
 * Run time: O(x^(2/3) / (log x)^2)
 * Memory usage: O(x^(1/2))
 */
int64_t primecount_nth_prime(int64_t n);

/*  Get the currently set number of threads */
int primecount_get_num_threads();

/*  Set the number of threads */
void primecount_set_num_threads(int num_threads);

/*
 * Largest number supported by primecount_pi_str(char* x).
 * @return 64-bit CPUs: max >= 10^27,
 *         32-bit CPUs: 2^63-1
 * The return type is a char* as primecount_get_max_x() may be a
 * 128-bit integer which is not supported by some compilers.
 */
const char* primecount_get_max_x();

/* Get the primecount version number, in the form “i.j” */
const char* primecount_version();

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
