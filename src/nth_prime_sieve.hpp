///
/// @file  nth_prime_sieve.hpp
/// @brief In the nth prime algorithm we first count the number of
///        primes up to an nth prime approximation. Next, we generate
///        primes using a special segmented sieve of Eratosthenes
///        algorithm with low memmory usage to find the actual nth
///        prime (which is close to the nth prime approximation).
///
///        Since we need to generate prime numbers close to the nth
///        prime which could potentially be as large as 10^30, we
///        cannot use the traditoinal segmented sieve of Eratosthenes
///        due to its O(n^(1/2)) memory usage. Therefore our
///        implementation uses a segment size of O(n^(1/3)) which
///        slightly deteriorates the runtime complexity of our
///        segmented sieve of Eratosthenes implementation. However,
///        our nth prime approximation is off by less than n^(1/2)
///        and therefore the slightly worse runtime complexity
///        of our sieving algorithm does not deteriorate the overall
///        runtime complexity of our nth prime algorithm.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef NTH_PRIME_SIEVE_HPP
#define NTH_PRIME_SIEVE_HPP

#include <int128_t.hpp>

namespace primecount {

int64_t nth_prime_sieve(int64_t n,
                        int64_t nth_prime_approx,
                        int64_t count_approx,
                        int threads);

#if defined(HAVE_INT128_T)

int128_t nth_prime_sieve(int128_t n,
                         int128_t nth_prime_approx,
                         int128_t count_approx,
                         int threads);

#endif

} // namespace

#endif
