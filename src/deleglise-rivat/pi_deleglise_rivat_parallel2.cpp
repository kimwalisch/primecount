///
/// @file  pi_deleglise_rivat_parallel2.cpp
/// @brief Parallel implementation of the Deleglise-Rivat prime
///        counting algorithm. Compared to
///        pi_deleglise_rivat_parallel1.cpp this version uses
///        compression (FactorTable & PiTable) to reduce the memory
///        usage.
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
#include <primecount-internal.hpp>
#include <aligned_vector.hpp>
#include <generate.hpp>
#include <pmath.hpp>
#include <PhiTiny.hpp>
#include <S1.hpp>

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
int64_t S2_trivial(int64_t x,
                   int64_t y,
                   int64_t z,
                   int64_t c,
                   PiTable& pi,
                   vector<int32_t>& primes,
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
  int64_t S2_total = 0;
  double time = get_wtime();

  // Find all trivial leaves: n = primes[b] * primes[l]
  // which satisfy phi(x / n), b - 1) = 1
  #pragma omp parallel for num_threads(threads) reduction(+: S2_total)
  for (int64_t b = max(c, pi_sqrtz + 1); b < pi_y; b++)
  {
    int64_t prime = primes[b];
    S2_total += pi_y - pi[max(x / (prime * prime), prime)];
  }

  if (print_status())
    print_result("S2_trivial", S2_total, time);

  return S2_total;
}

/// Calculate the contribution of the special leaves.
/// @pre y > 0 && c > 1
///
int64_t S2(int64_t x,
           int64_t y,
           int64_t z,
           int64_t c,
           vector<int32_t>& primes,
           FactorTable<uint16_t>& factors,
           int threads)
{
  int64_t S2_total = 0;
  int64_t limit = z + 1;
  threads = validate_threads(threads, limit);
  PiTable pi(y);

  S2_total += S2_trivial(x, y, z, c, pi, primes, threads);
  S2_total += S2_easy(x, y, z, c, pi, primes, threads);
  S2_total += S2_sieve(x, y, z, c, pi, primes, factors, threads);

  return S2_total;
}

/// alpha is a tuning factor which should grow like (log(x))^3
/// for the Deleglise-Rivat prime counting algorithm.
///
double compute_alpha(int64_t x)
{
  double d = (double) x;
  double alpha = log(d) * log(d) * log(d) / 1200;
  return in_between(1, alpha, iroot<6>(x));
}

} // namespace

namespace primecount {

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
int64_t pi_deleglise_rivat_parallel2(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  double alpha = compute_alpha(x);
  int64_t y = (int64_t) (alpha * iroot<3>(x));
  int64_t z = x / y;

  if (print_status())
  {
    cout << endl;
    cout << "=== pi_deleglise_rivat_parallel2(x) ===" << endl;
    cout << "pi(x) = S1 + S2 + pi(y) - 1 - P2" << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
    cout << "c = " << PhiTiny::max_a() << endl;
    cout << "threads = " << validate_threads(threads) << endl;
  }

  int64_t p2 = P2(x, y, threads);
  vector<int32_t> primes = generate_primes(y);
  FactorTable<uint16_t> factors(y);

  int64_t pi_y = pi_bsearch(primes, y);
  int64_t c = min(pi_y, PhiTiny::max_a());
  int64_t s1 = S1(x, y, c, primes[c], factors, threads);
  int64_t s2 = S2(x, y, z, c, primes, factors, threads);
  int64_t phi = s1 + s2;
  int64_t sum = phi + pi_y - 1 - p2;

  return sum;
}

} // namespace primecount
