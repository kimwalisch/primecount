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
///         1) Add a new option enum in CmdOptions.h.
///         2) Add your option to parseOptions() in CmdOptions.cpp.
///         3) Add your option to main() in main.cpp.
///         4) Document your option in help.cpp (--help option summary)
///            and in doc/primecount.txt (manpage).
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
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
  if (x > pstd::numeric_limits<int64_t>::max())
    throw primecount_error("x must be < 2^63");
  return (int64_t) x;
}

maxint_t AC(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  auto alpha = get_alpha_gourdon(x);
  double alpha_y = alpha.first;
  double alpha_z = alpha.second;
  maxint_t limit = get_max_x(alpha_y);

  if (x > limit)
    throw primecount_error("AC(x): x must be <= " + to_string(limit));

  int64_t x13 = iroot<3>(x);
  int64_t sqrtx = isqrt(x);
  int64_t y = (int64_t)(x13 * alpha_y);

  // x^(1/3) < y < x^(1/2)
  y = std::max(y, x13 + 1);
  y = std::min(y, sqrtx - 1);
  y = std::max(y, (int64_t) 1);

  int64_t k = PhiTiny::get_k(x);
  int64_t z = (int64_t)(y * alpha_z);

  // y <= z < x^(1/2)
  z = std::max(z, y);
  z = std::min(z, sqrtx - 1);
  z = std::max(z, (int64_t) 1);

  if (is_print())
    set_print_variables(true);

  if (x <= pstd::numeric_limits<int64_t>::max())
    return AC((int64_t) x, y, z, k, threads);
  else
    return AC(x, y, z, k, threads);
}

maxint_t B(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  auto alpha = get_alpha_gourdon(x);
  double alpha_y = alpha.first;
  maxint_t limit = get_max_x(alpha_y);

  if (x > limit)
    throw primecount_error("B(x): x must be <= " + to_string(limit));

  int64_t x13 = iroot<3>(x);
  int64_t sqrtx = isqrt(x);
  int64_t y = (int64_t)(x13 * alpha_y);

  // x^(1/3) < y < x^(1/2)
  y = std::max(y, x13 + 1);
  y = std::min(y, sqrtx - 1);
  y = std::max(y, (int64_t) 1);

  if (is_print())
    set_print_variables(true);

  if (x <= pstd::numeric_limits<int64_t>::max())
    return B((int64_t) x, y, threads);
  else
    return B(x, y, threads);
}

maxint_t D(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  auto alpha = get_alpha_gourdon(x);
  double alpha_y = alpha.first;
  double alpha_z = alpha.second;
  maxint_t limit = get_max_x(alpha_y);

  if (x > limit)
    throw primecount_error("D(x): x must be <= " + to_string(limit));

  int64_t x13 = iroot<3>(x);
  int64_t sqrtx = isqrt(x);
  int64_t y = (int64_t)(x13 * alpha_y);

  // x^(1/3) < y < x^(1/2)
  y = std::max(y, x13 + 1);
  y = std::min(y, sqrtx - 1);
  y = std::max(y, (int64_t) 1);

  int64_t k = PhiTiny::get_k(x);
  int64_t z = (int64_t)(y * alpha_z);

  // y <= z < x^(1/2)
  z = std::max(z, y);
  z = std::min(z, sqrtx - 1);
  z = std::max(z, (int64_t) 1);

  if (is_print())
    set_print_variables(true);

  if (x <= pstd::numeric_limits<int64_t>::max())
    return D((int64_t) x, y, z, k, (int64_t) Li(x), threads);
  else
    return D(x, y, z, k, Li(x), threads);
}

maxint_t Phi0(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  auto alpha = get_alpha_gourdon(x);
  double alpha_y = alpha.first;
  double alpha_z = alpha.second;
  maxint_t limit = get_max_x(alpha_y);

  if (x > limit)
    throw primecount_error("Phi0(x): x must be <= " + to_string(limit));

  int64_t x13 = iroot<3>(x);
  int64_t sqrtx = isqrt(x);
  int64_t y = (int64_t)(x13 * alpha_y);

  // x^(1/3) < y < x^(1/2)
  y = std::max(y, x13 + 1);
  y = std::min(y, sqrtx - 1);
  y = std::max(y, (int64_t) 1);

  int64_t k = PhiTiny::get_k(x);
  int64_t z = (int64_t)(y * alpha_z);

  // y <= z < x^(1/2)
  z = std::max(z, y);
  z = std::min(z, sqrtx - 1);
  z = std::max(z, (int64_t) 1);

  if (is_print())
    set_print_variables(true);

  if (x <= pstd::numeric_limits<int64_t>::max())
    return Phi0((int64_t) x, y, z, k, threads);
  else
    return Phi0(x, y, z, k, threads);
}

maxint_t Sigma(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  auto alpha = get_alpha_gourdon(x);
  double alpha_y = alpha.first;
  maxint_t limit = get_max_x(alpha_y);

  if (x > limit)
    throw primecount_error("Sigma(x): x must be <= " + to_string(limit));

  int64_t x13 = iroot<3>(x);
  int64_t sqrtx = isqrt(x);
  int64_t y = (int64_t)(x13 * alpha_y);

  // x^(1/3) < y < x^(1/2)
  y = std::max(y, x13 + 1);
  y = std::min(y, sqrtx - 1);
  y = std::max(y, (int64_t) 1);

  if (is_print())
    set_print_variables(true);

  if (x <= pstd::numeric_limits<int64_t>::max())
    return Sigma((int64_t) x, y, threads);
  else
    return Sigma(x, y, threads);
}

maxint_t P2(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  maxint_t limit = get_max_x(alpha);

  if (x > limit)
    throw primecount_error("P2(x): x must be <= " + to_string(limit));

  if (is_print())
    set_print_variables(true);

  int64_t y = (int64_t) (iroot<3>(x) * alpha);
  int64_t a = pi_noprint(y, threads);

  if (x <= pstd::numeric_limits<int64_t>::max())
    return P2((int64_t) x, y, a, threads);
  else
    return P2(x, y, a, threads);
}

maxint_t S1(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  maxint_t limit = get_max_x(alpha);

  if (x > limit)
    throw primecount_error("S1(x): x must be <= " + to_string(limit));

  if (is_print())
    set_print_variables(true);

  int64_t y = (int64_t) (iroot<3>(x) * alpha);
  int64_t c = PhiTiny::get_c(y);

  if (x <= pstd::numeric_limits<int64_t>::max())
    return S1((int64_t) x, y, c, threads);
  else
    return S1(x, y, c, threads);
}

maxint_t S2_trivial(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  maxint_t limit = get_max_x(alpha);

  if (x > limit)
    throw primecount_error("S2_trivial(x): x must be <= " + to_string(limit));

  if (is_print())
    set_print_variables(true);

  int64_t y = (int64_t) (iroot<3>(x) * alpha);
  int64_t z = (int64_t) (x / y);
  int64_t c = PhiTiny::get_c(y);

  if (x <= pstd::numeric_limits<int64_t>::max())
    return S2_trivial((int64_t) x, y, z, c, threads);
  else
    return S2_trivial(x, y, z, c, threads);
}

maxint_t S2_easy(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  maxint_t limit = get_max_x(alpha);

  if (x > limit)
    throw primecount_error("S2_easy(x): x must be <= " + to_string(limit));

  if (is_print())
    set_print_variables(true);

  int64_t y = (int64_t) (iroot<3>(x) * alpha);
  int64_t z = (int64_t) (x / y);
  int64_t c = PhiTiny::get_c(y);

  if (x <= pstd::numeric_limits<int64_t>::max())
    return S2_easy((int64_t) x, y, z, c, threads);
  else
    return S2_easy(x, y, z, c, threads);
}

maxint_t S2_hard(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  maxint_t limit = get_max_x(alpha);

  if (x > limit)
    throw primecount_error("S2_hard(x): x must be <= " + to_string(limit));

  if (is_print())
    set_print_variables(true);

  int64_t y = (int64_t) (iroot<3>(x) * alpha);
  int64_t z = (int64_t) (x / y);
  int64_t c = PhiTiny::get_c(y);

  if (x <= pstd::numeric_limits<int64_t>::max())
    return S2_hard((int64_t) x, y, z, c, (int64_t) Li(x), threads);
  else
    return S2_hard(x, y, z, c, Li(x), threads);
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
    auto threads = get_num_threads();
    maxint_t res = 0;

    switch (opts.option)
    {
      case OPTION_DEFAULT:
        res = pi(x, threads); break;
      case OPTION_DELEGLISE_RIVAT:
        res = pi_deleglise_rivat(x, threads); break;
      case OPTION_DELEGLISE_RIVAT_64:
        res = pi_deleglise_rivat_64(to_int64(x), threads); break;
      case OPTION_GOURDON:
        res = pi_gourdon(x, threads); break;
      case OPTION_GOURDON_64:
        res = pi_gourdon_64(to_int64(x), threads); break;
      case OPTION_LEGENDRE:
        res = pi_legendre(to_int64(x), threads); break;
      case OPTION_LEHMER:
        res = pi_lehmer(to_int64(x), threads); break;
      case OPTION_LMO:
        res = pi_lmo_parallel(to_int64(x), threads); break;
      case OPTION_LMO1:
        res = pi_lmo1(to_int64(x)); break;
      case OPTION_LMO2:
        res = pi_lmo2(to_int64(x)); break;
      case OPTION_LMO3:
        res = pi_lmo3(to_int64(x)); break;
      case OPTION_LMO4:
        res = pi_lmo4(to_int64(x)); break;
      case OPTION_LMO5:
        res = pi_lmo5(to_int64(x)); break;
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
      case OPTION_NTHPRIME:
        res = nth_prime(x, threads); break;
      case OPTION_NTHPRIME_64:
        res = nth_prime_64(x, threads); break;
      case OPTION_PHI:
        res = phi(to_int64(x), a, threads); break;
      case OPTION_P2:
        res = P2(x, threads); break;
      case OPTION_S1:
        res = S1(x, threads); break;
      case OPTION_S2_EASY:
        res = S2_easy(x, threads); break;
      case OPTION_S2_HARD:
        res = S2_hard(x, threads); break;
      case OPTION_S2_TRIVIAL:
        res = S2_trivial(x, threads); break;
      case OPTION_AC:
        res = AC(x, threads); break;
      case OPTION_B:
        res = B(x, threads); break;
      case OPTION_D:
        res = D(x, threads); break;
      case OPTION_PHI0:
        res = Phi0(x, threads); break;
      case OPTION_SIGMA:
        res = Sigma(x, threads); break;
#ifdef HAVE_INT128_T
      case OPTION_DELEGLISE_RIVAT_128:
        res = pi_deleglise_rivat_128(x, threads); break;
      case OPTION_GOURDON_128:
        res = pi_gourdon_128(x, threads); break;
      case OPTION_NTHPRIME_128:
        res = nth_prime_128(x, threads); break;
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
