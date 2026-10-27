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

#include <int128_t.hpp>
#include <isqrt.hpp>
#include <PiTable.hpp>
#include <Vector.hpp>

#include <stdint.h>

namespace primecount {

class Sieve;
class SemiprimeEventSorter;
struct SemiprimeEvent;

/// mertens_y[b] = sum_{n <= y, lpf(n) > primes[b]} mu(n)
/// for k < b <= pi[sqrt(y)]. Only entries with
/// primes[b]^2 <= y < primes[b]^3 are correct, these are
/// the entries used by the hyperbola method in D and AC.
///
template <typename Primes>
Vector<int64_t> generate_mertens_y(int64_t y,
                                   int64_t k,
                                   const Primes& primes,
                                   const PiTable& pi)
{
  int64_t pi_sqrty = pi[isqrt(y)];
  Vector<int64_t> mertens_y(pi_sqrty + 1);
  int64_t mertens = 1 - pi[y] + pi_sqrty;

  for (int64_t b = pi_sqrty; b > k; b--)
  {
    mertens_y[b] = mertens;
    mertens += pi[y / primes[b]] - b - 1;
  }

  return mertens_y;
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
