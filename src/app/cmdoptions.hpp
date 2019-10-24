///
/// @file  cmdoptions.hpp
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef CMDOPTIONS_HPP
#define CMDOPTIONS_HPP

#include <int128_t.hpp>
#include <stdint.h>

namespace primecount {

enum OptionID
{
  OPTION_ALPHA_Y,
  OPTION_ALPHA_Z,
  OPTION_DEFAULT,
  OPTION_GOURDON,
  OPTION_GOURDON_64,
  OPTION_GOURDON_128,
  OPTION_HELP,
  OPTION_LEGENDRE,
  OPTION_MEISSEL,
  OPTION_NTHPRIME,
  OPTION_NUMBER,
  OPTION_PRIMESIEVE,
  OPTION_LI,
  OPTION_LIINV,
  OPTION_RI,
  OPTION_RIINV,
  OPTION_PHI,
  OPTION_AC,
  OPTION_B,
  OPTION_D,
  OPTION_PHI0,
  OPTION_SIGMA,
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
  int option = OPTION_DEFAULT;
  bool time = false;
};

CmdOptions parseOptions(int, char**);

} // namespace

#endif
