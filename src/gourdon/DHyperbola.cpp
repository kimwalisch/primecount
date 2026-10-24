///
/// @file  DHyperbola.cpp
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

#include "DHyperbola.hpp"
#include "SemiprimeEventSorter.hpp"

#include <fast_div.hpp>
#include <int128_t.hpp>
#include <isqrt.hpp>
#include <macros.hpp>
#include <min.hpp>
#include <PiTable.hpp>
#include <sieve/Sieve.hpp>
#include <Vector.hpp>

#include <stdint.h>

namespace {

using namespace primecount;

/// Compute the special leaves x/(p*y) < m <= y of the prime
/// p = primes[b] in the current segment [low, high[ using the
/// hyperbola method, if sqrt(x/p) < y < p^3. Then max_m is lowered
/// to x/(p*y) and the caller processes the remaining leaves using
/// the factor table.
///
/// Let X = floor(x/p), L = xpy = floor(X/y), H = sqrt_xp =
/// floor(sqrt(X)) and Phi(t) = phi(t, b - 1). The leaves L < m <= y
/// contribute
///   M(H) Phi(H) - M(y) Phi(L)
///   + sum_{L < n <= H} a(n) E(X/n)
///   - sum_{L < n <= H} (a(n) + c(n)) Phi(X/n),
/// where a(n) = [lpf(n) >= p], c(n) = mu(n)[lpf(n) > p], E = Phi - M.
/// Here only the terms whose positions are >= p^2 are computed, see
/// easy_hyperbola_sum(). Prime n > p have a(n) + c(n) = 0. The
/// remaining terms, the events, have n = q*r with p <= q <= r and
/// weight a(n) + c(n) = 1 if q = p or q = r, and 2 otherwise.
/// Events are processed in decreasing n, so count() is called with
/// increasing positions X/n.
///
template <typename T, typename Primes>
ALWAYS_INLINE T hyperbola_sum(T xp,
                              int64_t p,
                              int64_t b,
                              int64_t y,
                              int64_t low,
                              int64_t high,
                              int64_t phi,
                              int64_t mertens_y,
                              const Primes& primes,
                              const PiTable& pi,
                              Sieve& sieve,
                              Vector<SemiprimeEvent>& events,
                              SemiprimeEventSorter& sorter,
                              int64_t& max_m)
{
  int64_t square = p * p;
  int64_t sqrt_xp = isqrt(xp);

  if (sqrt_xp >= y || square <= y / p)
    return 0;

  int64_t xpy = fast_div64(xp, y);
  int64_t query_low = max(low, square);
  max_m = min(max_m, xpy);
  T sum = 0;

  if (query_low <= xpy && xpy < high)
  {
    int64_t phi_xpy = phi + sieve.count(xpy - low);
    sum -= T(phi_xpy) * mertens_y;
  }

  if (query_low <= sqrt_xp && sqrt_xp < high)
  {
    int64_t phi_sqrt_xp = phi + sieve.count(sqrt_xp - low);
    sum += T(phi_sqrt_xp) * (phi_sqrt_xp - phi_minus_mertens(sqrt_xp, p, b, pi));
  }

  // n in (L, H] with position X/n in [query_low, high[
  int64_t min_n = min(fast_div(xp, high), sqrt_xp);
  int64_t max_n = min(fast_div64(xp, query_low), sqrt_xp);
  min_n = max(min_n, xpy);

  if (min_n < max_n)
  {
    sum += sum_phi_minus_mertens(xp, p, b, min_n, max_n, primes, pi);
    int64_t max_i = pi[isqrt(max_n)];
    events.clear();

    for (int64_t i = b; i <= max_i; i++)
    {
      int64_t q = primes[i];
      int64_t min_j = max(i, pi[min_n / q] + 1);
      int64_t max_j = pi[max_n / q];

      for (int64_t j = max_j; j >= min_j; j--)
      {
        int64_t n = q * primes[j];
        int64_t weight = (i != b && i != j) + 1;
        events.push_back(SemiprimeEvent{n, weight});
      }
    }

    sorter.sort(events, p, min_n, max_n);

    for (const SemiprimeEvent& event : events)
    {
      int64_t xpn = fast_div64(xp, event.n);
      int64_t phi_xpn = phi + sieve.count(xpn - low);
      sum += phi_minus_mertens(xpn, p, b, pi);
      sum -= T(event.weight) * phi_xpn;
    }
  }

  return sum;
}

} // namespace

namespace primecount {

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
                       int64_t& max_m)
{
  return ::hyperbola_sum(xp, p, b, y, low, high, phi, mertens_y,
                         primes, pi, sieve, events, sorter, max_m);
}

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
                       int64_t& max_m)
{
  return ::hyperbola_sum(xp, p, b, y, low, high, phi, mertens_y,
                         primes, pi, sieve, events, sorter, max_m);
}

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
                        int64_t& max_m)
{
  return ::hyperbola_sum(xp, p, b, y, low, high, phi, mertens_y,
                         primes, pi, sieve, events, sorter, max_m);
}

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
                        int64_t& max_m)
{
  return ::hyperbola_sum(xp, p, b, y, low, high, phi, mertens_y,
                         primes, pi, sieve, events, sorter, max_m);
}

#endif

} // namespace
