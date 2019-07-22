///
/// @file  cmdoptions.hpp
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef CMDOPTIONS_HPP
#define CMDOPTIONS_HPP

#include <int128_t.hpp>
#include <stdint.h>
#include <vector>

namespace primecount {

enum OptionID
{
  OPTION_ALPHA,
  OPTION_ALPHA_Y,
  OPTION_ALPHA_Z,
  OPTION_DELEGLISE_RIVAT,
  OPTION_DELEGLISE_RIVAT1,
  OPTION_DELEGLISE_RIVAT2,
  OPTION_DELEGLISE_RIVAT_PARALLEL1,
  OPTION_DELEGLISE_RIVAT_PARALLEL2,
  OPTION_GOURDON,
  OPTION_HELP,
  OPTION_LEGENDRE,
  OPTION_LEHMER,
  OPTION_LMO,
  OPTION_LMO1,
  OPTION_LMO2,
  OPTION_LMO3,
  OPTION_LMO4,
  OPTION_LMO5,
  OPTION_LMO_PARALLEL,
  OPTION_LI,
  OPTION_LIINV,
  OPTION_MEISSEL,
  OPTION_NTHPRIME,
  OPTION_NUMBER,
  OPTION_P2,
  OPTION_PHI,
  OPTION_PI,
  OPTION_PRIMESIEVE,
  OPTION_RI,
  OPTION_RIINV,
  OPTION_S1,
  OPTION_S2_EASY,
  OPTION_S2_HARD,
  OPTION_S2_TRIVIAL,
  OPTION_STATUS,
  OPTION_TEST,
  OPTION_TIME,
  OPTION_THREADS,
  OPTION_VERSION
};

struct CmdOptions
{
  maxint_t x = -1;
  int64_t a = -1;
  int option = OPTION_PI;
  bool time = false;
  std::vector<maxint_t> numbers;
};

CmdOptions parseOptions(int, char**);

} // namespace

#endif
