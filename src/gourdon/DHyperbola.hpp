///
/// @file  DHyperbola.hpp
/// @brief Compute the special leaves x/(p*m) with x/(p*y) < m <= y
///        of primecount's D formula using a weighted Dirichlet
///        hyperbola rearrangement, in which the sieve queries of the
///        prime leaves cancel.
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
#include <PiTable.hpp>
#include <Vector.hpp>

#include <stdint.h>

namespace primecount {

class Sieve;
class SemiprimeEventSorter;
struct SemiprimeEvent;

/// Precompute the boundary coefficients Phi(y) - E2(y) for the
/// eligible primes > min_prime, see DHyperbola.cpp.
///
Vector<int64_t> generate_hyperbola_coefficients(int64_t y,
                                                int64_t k,
                                                int64_t min_prime,
                                                const Vector<uint32_t>& primes,
                                                const PiTable& pi);

Vector<int64_t> generate_hyperbola_coefficients(int64_t y,
                                                int64_t k,
                                                int64_t min_prime,
                                                const Vector<int64_t>& primes,
                                                const PiTable& pi);

uint64_t hyperbola_sum(uint64_t xp,
                       int64_t p,
                       int64_t b,
                       int64_t y,
                       int64_t low,
                       int64_t high,
                       int64_t phi,
                       int64_t coefficient,
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
                        int64_t coefficient,
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
                        int64_t coefficient,
                        const Vector<int64_t>& primes,
                        const PiTable& pi,
                        Sieve& sieve,
                        Vector<SemiprimeEvent>& events,
                        SemiprimeEventSorter& sorter,
                        int64_t& max_m);

#endif

} // namespace

#endif
