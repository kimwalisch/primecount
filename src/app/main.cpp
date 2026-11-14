///
/// @file   main.cpp
/// @brief  Command-line option handling for the primecount
///         command-line application. The user's command-line options
///         are first parsed in CmdOptions.cpp and stored in a
///         CmdOptions object. Afterwards we execute the function
///         corresponding to the user's command-line options in the
///         main() function in main.cpp.
///
///         How to add a new command-line option:
///
///         1) Add a new option enum in CmdOptions.hpp.
///         2) Add your option to parseOptions() in CmdOptions.cpp.
///         3) Add your option to main() in main.cpp.
///         4) Document your option in help.cpp (--help option summary)
///            and in doc/primecount.txt (manpage).
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "CmdOptions.hpp"

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <gourdon.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <PhiTiny.hpp>
#include <PiTable.hpp>
#include <print.hpp>
#include <S.hpp>

#include <stdint.h>
#include <exception>
#include <iostream>
#include <string>

using namespace primecount;

namespace primecount {

int64_t to_int64(maxint_t x)
{
  if (x < pstd::numeric_limits<int64_t>::min() ||
      x > pstd::numeric_limits<int64_t>::max())
    throw primecount_error("the selected algorithm requires -2^63 <= x < 2^63");
  return (int64_t) x;
}

maxint_t AC(maxint_t x,
            int threads,
            double alpha)
{
  if (x < 1)
    return 0;

  alpha = get_alpha_gourdon(x, alpha);
  int64_t sqrtx = isqrt(x);
  int64_t x13 = iroot<3>(x);
  int64_t x14 = iroot<4>(x);
  int64_t y = (int64_t)(x13 * alpha);
  int64_t k = PhiTiny::get_k(x14);

  // x^(1/3) < y <= x^(1/2)
  y = std::max(y, x13 + 1);
  y = std::min(y, sqrtx);

#ifdef HAVE_INT128_T
  // x / y is stored in int64_t variables, hence x / y must be
  // <= 2^63-1. We use 2^62 as a safety buffer to protect
  // against overflows in calculations derived from x / y.
  if (x > (int128_t(y) << 62))
    throw primecount_error("AC(x): x is too large");
#endif

  if (is_print())
    set_print_variables(true);

  const PiTable pi(y, threads);

  if (x <= pstd::numeric_limits<int64_t>::max())
    return AC((int64_t) x, y, k, pi, threads);
  else
    return AC(x, y, k, pi, threads);
}

maxint_t B(maxint_t x,
           int threads,
           double alpha)
{
  if (x < 1)
    return 0;

  alpha = get_alpha_gourdon(x, alpha);
  int64_t x13 = iroot<3>(x);
  int64_t sqrtx = isqrt(x);
  int64_t y = (int64_t)(x13 * alpha);

  // x^(1/3) < y <= x^(1/2)
  y = std::max(y, x13 + 1);
  y = std::min(y, sqrtx);

#ifdef HAVE_INT128_T
  // x / y is stored in int64_t variables, hence x / y must be
  // <= 2^63-1. We use 2^62 as a safety buffer to protect
  // against overflows in calculations derived from x / y.
  if (x > (int128_t(y) << 62))
    throw primecount_error("B(x): x is too large");
#endif

  if (is_print())
    set_print_variables(true);

  if (x <= pstd::numeric_limits<int64_t>::max())
    return B((int64_t) x, y, threads);
  else
    return B(x, y, threads);
}

maxint_t D(maxint_t x,
           int threads,
           double alpha)
{
  if (x < 1)
    return 0;

  alpha = get_alpha_gourdon(x, alpha);
  int64_t sqrtx = isqrt(x);
  int64_t x13 = iroot<3>(x);
  int64_t x14 = iroot<4>(x);
  int64_t y = (int64_t)(x13 * alpha);
  int64_t k = PhiTiny::get_k(x14);

  // x^(1/3) < y <= x^(1/2)
  y = std::max(y, x13 + 1);
  y = std::min(y, sqrtx);

#ifdef HAVE_INT128_T
  // x / y is stored in int64_t variables, hence x / y must be
  // <= 2^63-1. We use 2^62 as a safety buffer to protect
  // against overflows in calculations derived from x / y.
  if (x > (int128_t(y) << 62))
    throw primecount_error("D(x): x is too large");
#endif

  if (is_print())
    set_print_variables(true);

  const PiTable pi(y, threads);

  if (x <= pstd::numeric_limits<int64_t>::max())
    return D((int64_t) x, y, k, pi, threads);
  else
    return D(x, y, k, pi, threads);
}

maxint_t Phi0(maxint_t x,
              int threads,
              double alpha)
{
  if (x < 1)
    return 0;

  alpha = get_alpha_gourdon(x, alpha);
  int64_t sqrtx = isqrt(x);
  int64_t x13 = iroot<3>(x);
  int64_t x14 = iroot<4>(x);
  int64_t y = (int64_t)(x13 * alpha);
  int64_t k = PhiTiny::get_k(x14);

  // x^(1/3) < y <= x^(1/2)
  y = std::max(y, x13 + 1);
  y = std::min(y, sqrtx);

#ifdef HAVE_INT128_T
  // x / y is stored in int64_t variables, hence x / y must be
  // <= 2^63-1. We use 2^62 as a safety buffer to protect
  // against overflows in calculations derived from x / y.
  if (x > (int128_t(y) << 62))
    throw primecount_error("Phi0(x): x is too large");
#endif

  if (is_print())
    set_print_variables(true);

  const PiTable pi(y, threads);

  if (x <= pstd::numeric_limits<int64_t>::max())
    return Phi0((int64_t) x, y, k, pi, threads);
  else
    return Phi0(x, y, k, pi, threads);
}

maxint_t Sigma(maxint_t x,
               int threads,
               double alpha)
{
  if (x < 1)
    return 0;

  alpha = get_alpha_gourdon(x, alpha);
  int64_t x13 = iroot<3>(x);
  int64_t sqrtx = isqrt(x);
  int64_t y = (int64_t)(x13 * alpha);

  // x^(1/3) < y <= x^(1/2)
  y = std::max(y, x13 + 1);
  y = std::min(y, sqrtx);

#ifdef HAVE_INT128_T
  // x / y is stored in int64_t variables, hence x / y must be
  // <= 2^63-1. We use 2^62 as a safety buffer to protect
  // against overflows in calculations derived from x / y.
  if (x > (int128_t(y) << 62))
    throw primecount_error("Sigma(x): x is too large");
#endif

  if (is_print())
    set_print_variables(true);

  const PiTable pi(y, threads);

  if (x <= pstd::numeric_limits<int64_t>::max())
    return Sigma((int64_t) x, y, pi, threads);
  else
    return Sigma(x, y, pi, threads);
}

maxint_t P2(maxint_t x,
            int threads,
            double alpha)
{
  if (x < 1)
    return 0;

  alpha = get_alpha_deleglise_rivat(x, alpha);
  int64_t y = (int64_t) (iroot<3>(x) * alpha);

#ifdef HAVE_INT128_T
  // x / y is stored in int64_t variables, hence x / y must be
  // <= 2^63-1. We use 2^62 as a safety buffer to protect
  // against overflows in calculations derived from x / y.
  if (x > (int128_t(y) << 62))
    throw primecount_error("P2(x): x is too large");
#endif

  if (is_print())
    set_print_variables(true);

  int64_t a = pi_noprint(y, threads);

  if (x <= pstd::numeric_limits<int64_t>::max())
    return P2((int64_t) x, y, a, threads);
  else
    return P2(x, y, a, threads);
}

maxint_t S1(maxint_t x,
            int threads,
            double alpha)
{
  if (x < 1)
    return 0;

  alpha = get_alpha_deleglise_rivat(x, alpha);
  int64_t y = (int64_t) (iroot<3>(x) * alpha);

#ifdef HAVE_INT128_T
  // x / y is stored in int64_t variables, hence x / y must be
  // <= 2^63-1. We use 2^62 as a safety buffer to protect
  // against overflows in calculations derived from x / y.
  if (x > (int128_t(y) << 62))
    throw primecount_error("S1(x): x is too large");
#endif

  if (is_print())
    set_print_variables(true);

  int64_t k = PhiTiny::get_k(y);

  if (x <= pstd::numeric_limits<int64_t>::max())
    return S1((int64_t) x, y, k, threads);
  else
    return S1(x, y, k, threads);
}

maxint_t S2_trivial(maxint_t x,
                    int threads,
                    double alpha)
{
  if (x < 1)
    return 0;

  alpha = get_alpha_deleglise_rivat(x, alpha);
  int64_t y = (int64_t) (iroot<3>(x) * alpha);

#ifdef HAVE_INT128_T
  // x / y is stored in int64_t variables, hence x / y must be
  // <= 2^63-1. We use 2^62 as a safety buffer to protect
  // against overflows in calculations derived from x / y.
  if (x > (int128_t(y) << 62))
    throw primecount_error("S2_trivial(x): x is too large");
#endif

  int64_t xy = (int64_t) (x / y);
  int64_t k = PhiTiny::get_k(y);

  if (is_print())
    set_print_variables(true);

  if (x <= pstd::numeric_limits<int64_t>::max())
    return S2_trivial((int64_t) x, y, xy, k, threads);
  else
    return S2_trivial(x, y, xy, k, threads);
}

maxint_t S2_easy(maxint_t x,
                 int threads,
                 double alpha)
{
  if (x < 1)
    return 0;

  alpha = get_alpha_deleglise_rivat(x, alpha);
  int64_t y = (int64_t) (iroot<3>(x) * alpha);

#ifdef HAVE_INT128_T
  // x / y is stored in int64_t variables, hence x / y must be
  // <= 2^63-1. We use 2^62 as a safety buffer to protect
  // against overflows in calculations derived from x / y.
  if (x > (int128_t(y) << 62))
    throw primecount_error("S2_easy(x): x is too large");
#endif

  int64_t xy = (int64_t) (x / y);
  int64_t k = PhiTiny::get_k(y);

  if (is_print())
    set_print_variables(true);

  if (x <= pstd::numeric_limits<int64_t>::max())
    return S2_easy((int64_t) x, y, xy, k, threads);
  else
    return S2_easy(x, y, xy, k, threads);
}

maxint_t S2_hard(maxint_t x,
                 int threads,
                 double alpha)
{
  if (x < 1)
    return 0;

  alpha = get_alpha_deleglise_rivat(x, alpha);
  int64_t y = (int64_t) (iroot<3>(x) * alpha);

#ifdef HAVE_INT128_T
  // x / y is stored in int64_t variables, hence x / y must be
  // <= 2^63-1. We use 2^62 as a safety buffer to protect
  // against overflows in calculations derived from x / y.
  if (x > (int128_t(y) << 62))
    throw primecount_error("S2_hard(x): x is too large");
#endif

  int64_t xy = (int64_t) (x / y);
  int64_t k = PhiTiny::get_k(y);

  if (is_print())
    set_print_variables(true);

  if (x <= pstd::numeric_limits<int64_t>::max())
    return S2_hard((int64_t) x, y, xy, k, threads);
  else
    return S2_hard(x, y, xy, k, threads);
}

} // namespace

int main (int argc, char* argv[])
{
  try
  {
    CmdOptions opts = parseOptions(argc, argv);
    double time = get_time();

    auto x = opts.x;
    auto a = opts.a;
    auto alpha = opts.alpha;
    auto threads = get_num_threads();
    maxint_t res = 0;

    switch (opts.option)
    {
      case OPTION_DEFAULT:
        res = pi(x, threads, alpha); break;
      case OPTION_DELEGLISE_RIVAT:
        res = pi_deleglise_rivat(x, threads, alpha); break;
      case OPTION_DELEGLISE_RIVAT_64:
        res = pi_deleglise_rivat_64(to_int64(x), threads, alpha); break;
      case OPTION_GOURDON:
        res = pi_gourdon(x, threads, alpha); break;
      case OPTION_GOURDON_64:
        res = pi_gourdon_64(to_int64(x), threads, alpha); break;
      case OPTION_LEGENDRE:
        res = pi_legendre(to_int64(x), threads); break;
      case OPTION_LEHMER:
        res = pi_lehmer(to_int64(x), threads); break;
      case OPTION_LMO:
        res = pi_lmo_parallel(to_int64(x), threads, alpha); break;
      case OPTION_LMO1:
        res = pi_lmo1(to_int64(x)); break;
      case OPTION_LMO2:
        res = pi_lmo2(to_int64(x), alpha); break;
      case OPTION_LMO3:
        res = pi_lmo3(to_int64(x), alpha); break;
      case OPTION_LMO4:
        res = pi_lmo4(to_int64(x), alpha); break;
      case OPTION_LMO5:
        res = pi_lmo5(to_int64(x), alpha); break;
      case OPTION_MEISSEL:
        res = pi_meissel(to_int64(x), threads); break;
      case OPTION_PRIMESIEVE:
        res = pi_primesieve(to_int64(x)); break;
      case OPTION_LI:
        res = Li(x); break;
      case OPTION_LIINV:
        res = Li_inverse(x); break;
      case OPTION_R:
        res = RiemannR(x); break;
      case OPTION_R_INVERSE:
        res = RiemannR_inverse(x); break;
      case OPTION_R_PSI:
        res = RiemannR_psi(x); break;
      case OPTION_R_PSI_INVERSE:
        res = RiemannR_psi_inverse(x); break;
      case OPTION_NTHPRIME:
        res = nth_prime(x, threads, alpha); break;
      case OPTION_NTHPRIME_64:
        res = nth_prime_64(to_int64(x), threads, alpha); break;
      case OPTION_PHI:
        res = phi(to_int64(x), a, threads); break;
      case OPTION_P2:
        res = P2(x, threads, alpha); break;
      case OPTION_S1:
        res = S1(x, threads, alpha); break;
      case OPTION_S2_EASY:
        res = S2_easy(x, threads, alpha); break;
      case OPTION_S2_HARD:
        res = S2_hard(x, threads, alpha); break;
      case OPTION_S2_TRIVIAL:
        res = S2_trivial(x, threads, alpha); break;
      case OPTION_AC:
        res = AC(x, threads, alpha); break;
      case OPTION_B:
        res = B(x, threads, alpha); break;
      case OPTION_D:
        res = D(x, threads, alpha); break;
      case OPTION_PHI0:
        res = Phi0(x, threads, alpha); break;
      case OPTION_SIGMA:
        res = Sigma(x, threads, alpha); break;
#ifdef HAVE_INT128_T
      case OPTION_DELEGLISE_RIVAT_128:
        res = pi_deleglise_rivat_128(x, threads, alpha); break;
      case OPTION_GOURDON_128:
        res = pi_gourdon_128(x, threads, alpha); break;
      case OPTION_NTHPRIME_128:
        res = nth_prime_128(x, threads, alpha); break;
#endif
    }

    if (is_print_combined_result())
    {
      // Add empty line after last partial formula
      if (is_print())
        std::cout << std::endl;

      std::cout << res << std::endl;

      if (opts.time)
        print_seconds(get_time() - time);
    }
  }
  catch (std::exception& e)
  {
    std::cerr << "primecount: " << e.what() << std::endl
              << "Try 'primecount --help' for more information." << std::endl;
    return 1;
  }

  return 0;
}
