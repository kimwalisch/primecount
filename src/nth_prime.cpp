///
/// @file  nth_prime.cpp
/// @brief Find the nth prime.
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "nth_prime_sieve.hpp"
#include "PiTable.hpp"

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <Vector.hpp>
#include <imath.hpp>
#include <macros.hpp>

#include <stdint.h>
#include <string>

namespace {

using namespace primecount;

// Number of primes < 2^63
constexpr int64_t max_n_int64 = 216289611853439384ll;

// primes[1] = 2, primes[2] = 3, ...
const Array<int16_t, 170> primes =
{
    0,   2,   3,   5,   7,  11,  13,  17,  19,  23,
   29,  31,  37,  41,  43,  47,  53,  59,  61,  67,
   71,  73,  79,  83,  89,  97, 101, 103, 107, 109,
  113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
  173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
  229, 233, 239, 241, 251, 257, 263, 269, 271, 277,
  281, 283, 293, 307, 311, 313, 317, 331, 337, 347,
  349, 353, 359, 367, 373, 379, 383, 389, 397, 401,
  409, 419, 421, 431, 433, 439, 443, 449, 457, 461,
  463, 467, 479, 487, 491, 499, 503, 509, 521, 523,
  541, 547, 557, 563, 569, 571, 577, 587, 593, 599,
  601, 607, 613, 617, 619, 631, 641, 643, 647, 653,
  659, 661, 673, 677, 683, 691, 701, 709, 719, 727,
  733, 739, 743, 751, 757, 761, 769, 773, 787, 797,
  809, 811, 821, 823, 827, 829, 839, 853, 857, 859,
  863, 877, 881, 883, 887, 907, 911, 919, 929, 937,
  941, 947, 953, 967, 971, 977, 983, 991, 997, 1009
};

/// Find the nth prime using binary search
/// and a PrimePi(x) lookup table.
/// Run time: O(log2(n))
///
int64_t binary_search_nth_prime(int64_t n)
{
  int64_t low = n * 2;
  int64_t hi = PiTable::max_cached();

  ASSERT(low < hi);
  ASSERT(n >= 1);
  ASSERT(n >= PiTable::pi_cache(low));
  ASSERT(n <= PiTable::pi_cache(hi));

  while (low < hi)
  {
    int64_t mid = low + (hi - low) / 2;
    if (PiTable::pi_cache(mid) < n)
      low = mid + 1;
    else
      hi = mid;
  }

  return low;
}

} // namespace

namespace primecount {

/// Find the nth prime using the prime counting function
/// and the segmented sieve of Eratosthenes.
/// Run time: O(x^(2/3) / (log x)^2)
/// Memory usage: O(x^(1/3) * (log x)^3)
///
int64_t nth_prime_64(int64_t n, int threads)
{
  if_unlikely(n < 1)
    throw primecount_error("nth_prime(n): n must be >= 1");
  if_unlikely(n > max_n_int64)
    throw primecount_error("nth_prime(n): n must be <= " + std::to_string(max_n_int64));

  // For tiny n <= 169
  if (n < (int64_t) primes.size())
    return primes[n];

  // For small n <= 3314
  if (n <= PiTable::pi_cache(PiTable::max_cached()))
    return binary_search_nth_prime(n);

  // Closely approximate the nth prime using the inverse
  // Riemann R function and then count the primes up to this
  // approximation using the prime counting function.
  int64_t prime_approx = RiemannR_inverse(n);
  int64_t count_approx = pi(prime_approx, threads);

  constexpr bool forward = true;
  constexpr bool backward = false;

  // Use multi-threaded NthPrimeSieve for
  // large nth prime computations.
  if (threads > 1 &&
      prime_approx > (int64_t) 1e16)
  {
    // Here we are very close to the nth prime < sqrt(nth_prime),
    // we use a prime sieve to find the actual nth prime.
    if (count_approx < n)
      return nth_prime_sieve<forward>(n - count_approx, prime_approx + 1, threads);
    else
      return nth_prime_sieve<backward>(1 + count_approx - n, prime_approx, threads);
  }

  int64_t avg_prime_gap = ilog(prime_approx) + 2;
  int64_t prime = -1;

  // Here we are very close to the nth prime < sqrt(nth_prime),
  // we simply iterate over the primes until we find it.
  if (count_approx < n)
  {
    uint64_t start = prime_approx + 1;
    uint64_t stop = start + (n - count_approx) * avg_prime_gap;
    primesieve::iterator iter(start, stop);
    for (int64_t i = count_approx; i < n; i++)
      prime = iter.next_prime();
  }
  else // if (count_approx >= n)
  {
    uint64_t start = prime_approx;
    uint64_t stop = start - (count_approx - n) * avg_prime_gap;
    primesieve::iterator iter(start, stop);
    for (int64_t i = count_approx; i >= n; i--)
      prime = iter.prev_prime();
  }

  return prime;
}

#if defined(HAVE_INT128_T)

/// Find the nth prime using the prime counting function
/// and the segmented sieve of Eratosthenes.
/// Run time: O(x^(2/3) / (log x)^2)
/// Memory usage: O(x^(1/3) * (log x)^3)
///
int128_t nth_prime_128(int128_t n, int threads)
{
  if_unlikely(n < 1)
    throw primecount_error("nth_prime(n): n must be >= 1");

  // For tiny n <= 169
  if (n < (int64_t) primes.size())
    return primes[n];

  // For small n <= 3314
  if (n <= PiTable::pi_cache(PiTable::max_cached()))
    return binary_search_nth_prime(n);

  // Closely approximate the nth prime using the inverse
  // Riemann R function and then count the primes up to this
  // approximation using the prime counting function.
  int128_t prime_approx = RiemannR_inverse(n);
  int128_t count_approx = pi(prime_approx, threads);

  constexpr bool forward = true;
  constexpr bool backward = false;

  // Here we are very close to the nth prime < sqrt(nth_prime),
  // we use a prime sieve to find the actual nth prime.
  if (count_approx < n)
    return nth_prime_sieve<forward>(n - count_approx, prime_approx + 1, threads);
  else
    return nth_prime_sieve<backward>(1 + count_approx - n, prime_approx, threads);
}

#endif

} // namespace
