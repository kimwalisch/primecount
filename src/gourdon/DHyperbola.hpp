///
/// @file  DHyperbola.hpp
/// @brief Compute the special leaves x/(p*m) with x/(p*y) < m <= y
///        of Gourdon's D formula using a weighted Dirichlet hyperbola
///        rearrangement, in which the sieve queries of the prime
///        leaves cancel.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
/// Copyright (C) 2026 Chopper Lin, <choppers@126.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef D_HYPERBOLA_HPP
#define D_HYPERBOLA_HPP

#include <fast_div.hpp>
#include <int128_t.hpp>
#include <isqrt.hpp>
#include <macros.hpp>
#include <PiTable.hpp>
#include <Vector.hpp>

#include <stdint.h>

namespace primecount {

class Sieve;
class SemiprimeEventSorter;

struct SemiprimeEvent
{
  int64_t n;
  int64_t weight;
};

/// Returns phi(t, b - 1) - M(t) for p = primes[b] and p^2 <= t < p^3,
/// where M(t) = sum_{n <= t, lpf(n) > p} mu(n).
///
ALWAYS_INLINE int64_t phi_minus_mertens(int64_t t,
                                        int64_t p,
                                        int64_t b,
                                        const PiTable& pi)
{
  return 2 * pi[t] + pi[t / p] + pi[isqrt(t)] - 4 * b + 2;
}

/// Sum pi(xp / primes[j]) for i <= j <= last
template <typename T, typename Primes>
ALWAYS_INLINE T sum_pi(T xp,
                       int64_t i,
                       int64_t last,
                       const Primes& primes,
                       const PiTable& pi)
{
  T sum = 0;

  // Unroll loop to increase instruction level parallelism
  for (; i + 3 <= last; i += 4)
  {
    int64_t xpq0 = fast_div64(xp, primes[i]);
    int64_t xpq1 = fast_div64(xp, primes[i + 1]);
    int64_t xpq2 = fast_div64(xp, primes[i + 2]);
    int64_t xpq3 = fast_div64(xp, primes[i + 3]);
    sum += pi[xpq0] + pi[xpq1] + pi[xpq2] + pi[xpq3];
  }

  NO_UNROLL_LOOP
  for (; i <= last; i++)
    sum += pi[fast_div64(xp, primes[i])];

  return sum;
}

/// Sum phi_minus_mertens(X/q) over the primes q in (min_n, max_n],
/// where X = x/p and all X/q are in [p^2, p^3[. Gourdon's inversion
/// equality turns the sums of pi(X/(p*q)) and pi(sqrt(X/q)) into sums
/// over the reflected primes r, avoiding two pi lookups per prime q.
///
template <typename T, typename Primes>
T sum_phi_minus_mertens(T xp,
                        int64_t p,
                        int64_t b,
                        int64_t min_n,
                        int64_t max_n,
                        const Primes& primes,
                        const PiTable& pi)
{
  int64_t pi_min_n = pi[min_n];
  int64_t pi_max_n = pi[max_n];

  if (pi_min_n == pi_max_n)
    return 0;

  T xpp = xp / p;
  int64_t min_r = fast_div64(xpp, max_n);
  int64_t max_r = fast_div64(xpp, min_n + 1);
  int64_t min_r2 = isqrt(fast_div64(xp, max_n));
  int64_t max_r2 = isqrt(fast_div64(xp, min_n + 1));

  T sum = sum_pi(xp, pi_min_n + 1, pi_max_n, primes, pi) * 2;
  sum += sum_pi(xpp, pi[min_r] + 1, pi[max_r], primes, pi);
  sum += T(pi_max_n) * pi[min_r] - T(pi_min_n) * pi[max_r];
  sum += T(pi_max_n) * pi[min_r2] - T(pi_min_n) * pi[max_r2];
  sum -= T(pi_max_n - pi_min_n) * (4 * b - 2);

  for (int64_t j = pi[min_r2] + 1, last = pi[max_r2]; j <= last; j++)
  {
    int64_t r = primes[j];
    sum += pi[fast_div64(xp, r * r)];
  }

  return sum;
}

/// Compute the hyperbola terms whose positions are < p^2. These only
/// need PiTable lookups, and D_thread() skips the prime p in threads
/// whose sieve limit is < p^2, hence D.cpp computes them once per
/// prime and hyperbola_sum() skips them. Only for the primes
/// (x/y)^(1/3) < p <= sqrt(y): x < p^3 * y <= p * y^2 and y^2 <= x
/// imply sqrt(x/p) < y < p^3, and their leaves include easy leaves.
///
template <typename T, typename Primes>
T easy_hyperbola_sum(T xp,
                     int64_t p,
                     int64_t b,
                     int64_t y,
                     int64_t mertens_y,
                     const Primes& primes,
                     const PiTable& pi)
{
  int64_t square = p * p;
  int64_t x_div_p3 = fast_div64(xp, square);
  int64_t xpy = fast_div64(xp, y);
  int64_t sqrt_xp = isqrt(xp);
  T sum = 0;

  // xpy = X/y < p^2, phi(xpy, b - 1) counts
  // 1 and the primes in [p, xpy].
  int64_t phi_xpy = pi[xpy] - b + 2;
  sum -= T(mertens_y) * phi_xpy;

  // If H >= p^2, only the L term has a position < p^2
  if (sqrt_xp < square)
  {
    int64_t phi_sqrt_xp = pi[sqrt_xp] - b + 2;
    int64_t mertens_sqrt_xp = 1 - (pi[sqrt_xp] - b);
    sum += T(mertens_sqrt_xp) * phi_sqrt_xp;

    // Primes n in (X/p^2, H] have X/n < p^2 (X/p^2 >= L as p^2 <= y):
    // E(X/n) = phi(X/n, b - 1) - M(X/n) = 2 * pi[X/n] - 2b + 1.
    int64_t pi_min_n = pi[x_div_p3];
    int64_t pi_max_n = pi[sqrt_xp];
    sum += sum_pi(xp, pi_min_n + 1, pi_max_n, primes, pi) * 2;
    sum -= T(pi_max_n - pi_min_n) * (2 * b - 1);
  }

  return sum;
}

uint64_t hyperbola_sum(uint64_t xp,
                       int64_t p,
                       int64_t b,
                       int64_t y,
                       int64_t low,
                       int64_t high,
                       int64_t phi,
                       int64_t mertens_y,
                       const Vector<uint32_t>& primes,
                       const PiTable& pi,
                       Sieve& sieve,
                       Vector<SemiprimeEvent>& events,
                       SemiprimeEventSorter& sorter,
                       int64_t& max_m);

uint64_t hyperbola_sum(uint64_t xp,
                       int64_t p,
                       int64_t b,
                       int64_t y,
                       int64_t low,
                       int64_t high,
                       int64_t phi,
                       int64_t mertens_y,
                       const Vector<int64_t>& primes,
                       const PiTable& pi,
                       Sieve& sieve,
                       Vector<SemiprimeEvent>& events,
                       SemiprimeEventSorter& sorter,
                       int64_t& max_m);

#if defined(HAVE_INT128_T)

uint128_t hyperbola_sum(uint128_t xp,
                        int64_t p,
                        int64_t b,
                        int64_t y,
                        int64_t low,
                        int64_t high,
                        int64_t phi,
                        int64_t mertens_y,
                        const Vector<uint32_t>& primes,
                        const PiTable& pi,
                        Sieve& sieve,
                        Vector<SemiprimeEvent>& events,
                        SemiprimeEventSorter& sorter,
                        int64_t& max_m);

uint128_t hyperbola_sum(uint128_t xp,
                        int64_t p,
                        int64_t b,
                        int64_t y,
                        int64_t low,
                        int64_t high,
                        int64_t phi,
                        int64_t mertens_y,
                        const Vector<int64_t>& primes,
                        const PiTable& pi,
                        Sieve& sieve,
                        Vector<SemiprimeEvent>& events,
                        SemiprimeEventSorter& sorter,
                        int64_t& max_m);

#endif

} // namespace

#endif
