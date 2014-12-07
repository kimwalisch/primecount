///
/// @file  pi_deleglise_rivat_parallel3.cpp
/// @brief Parallel implementation of the Deleglise-Rivat prime
///        counting algorithm. This implementation is identical to
///        pi_deleglise_rivat_parallel2(x) but uses 128-bit integers.
///
///        This implementation is based on the paper:
///        Tomás Oliveira e Silva, Computing pi(x): the combinatorial
///        method, Revista do DETUA, vol. 4, no. 6, March 2006,
///        pp. 759-768.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <FactorTable.hpp>
#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <generate.hpp>
#include <pmath.hpp>
#include <PhiTiny.hpp>
#include <int128.hpp>
#include <S1.hpp>
#include "S2.hpp"

#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;

namespace {

/// Calculate the contribution of the trivial leaves.
///
template <typename P>
int128_t S2_trivial(uint128_t x,
                   int64_t y,
                   int64_t z,
                   int64_t c,
                   PiTable& pi,
                   vector<P>& primes,
                   int threads)
{
  if (print_status())
  {
    cout << endl;
    cout << "=== S2_trivial(x, y) ===" << endl;
    cout << "Computation of the trivial special leaves" << endl;
  }

  int64_t pi_y = pi[y];
  int64_t pi_sqrtz = pi[min(isqrt(z), y)];
  int128_t S2_total = 0;
  double time = get_wtime();

  // Find all trivial leaves: n = primes[b] * primes[l]
  // which satisfy phi(x / n), b - 1) = 1
  #pragma omp parallel for num_threads(threads) reduction(+: S2_total)
  for (int64_t b = max(c, pi_sqrtz + 1); b < pi_y; b++)
  {
    uint128_t prime = primes[b];
    uint64_t xn = (uint64_t) max(x / (prime * prime), prime);
    S2_total += pi_y - pi[xn];
  }

  if (print_status())
    print_result("S2_trivial", S2_total, time);

  return S2_total;
}

/// Calculate the contribution of the special leaves.
/// @pre y > 0 && c > 1
///
template <typename P, typename F>
int128_t S2(uint128_t x,
            int64_t y,
            int64_t z,
            int64_t c,
            vector<P>& primes,
            FactorTable<F>& factors,
            int threads)
{
  int128_t S2_total = 0;

  threads = validate_threads(threads, z);
  PiTable pi(y);

  S2_total += S2_trivial(x, y, z, c, pi, primes, threads);
  S2_total += S2_easy(x, y, z, c, pi, primes, threads);
  S2_total += S2_sieve(x, y, z, c, pi, primes, factors, threads);

  return S2_total;
}

/// alpha is a tuning factor which should grow like (log(x))^3
/// for the Deleglise-Rivat prime counting algorithm.
///
double compute_alpha(int128_t x)
{
  double d = (double) x;
  double alpha = log(d) * log(d) * log(d) / 1000;
  return in_between(1, alpha, iroot<6>(x));
}

} // namespace

namespace primecount {

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
int128_t pi_deleglise_rivat_parallel3(int128_t x, int threads)
{
  if (x < 2)
    return 0;

  if (x > to_maxint(primecount::max()))
    throw primecount_error("pi(x): x must be <= " + max());

  double alpha = compute_alpha(x);
  int64_t y = (int64_t) (alpha * iroot<3>(x));
  int64_t z = (int64_t) (x / y);
  int64_t pi_y;
  int64_t c;

  if (print_status())
  {
    cout << endl;
    cout << "=== pi_deleglise_rivat_parallel3(x) ===" << endl;
    cout << "pi(x) = S1 + S2 + pi(y) - 1 - P2" << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
    cout << "c = " << PhiTiny::max_a() << endl;
    cout << "threads = " << validate_threads(threads) << endl;
  }

  int128_t p2 = P2(x, y, threads);
  int128_t s1, s2;

  if (y <= FactorTable<uint16_t>::max())
  {
    // if y < 2^32 we can use 32-bit primes and a 16-bit FactorTable
    // which uses ~ (y / 2) bytes of memory

    vector<uint32_t> primes = generate_primes<uint32_t>(y);
    FactorTable<uint16_t> factors(y);

    pi_y = primes.size() - 1;
    c = min(pi_y, PhiTiny::max_a());
    s1 = S1(x, y, c, primes[c], factors, threads);
    s2 = S2(x, y, z, c, primes, factors, threads);
  }
  else
  {
    // if y >= 2^32 we need to use 64-bit primes and a 32-bit
    // FactorTable which uses ~ y bytes of memory

    vector<int64_t> primes = generate_primes<int64_t>(y);
    FactorTable<uint32_t> factors(y);

    pi_y = primes.size() - 1;
    c = min(pi_y, PhiTiny::max_a());
    s1 = S1(x, y, c, primes[c], factors, threads);
    s2 = S2(x, y, z, c, primes, factors, threads);
  }

  int128_t phi = s1 + s2;
  int128_t sum = phi + pi_y - 1 - p2;

  return sum;
}

} // namespace
