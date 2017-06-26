///
/// @file  cmdoptions.hpp
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef CMDOPTIONS_HPP
#define CMDOPTIONS_HPP

#include <primecount.hpp>
#include <int128_t.hpp>
#include <stdint.h>
#include <string>

namespace primecount {

enum OptionID
{
  OPTION_ALPHA,
  OPTION_BACKUP,
  OPTION_DELEGLISE_RIVAT,
  OPTION_DELEGLISE_RIVAT1,
  OPTION_DELEGLISE_RIVAT2,
  OPTION_DELEGLISE_RIVAT_PARALLEL1,
  OPTION_DELEGLISE_RIVAT_PARALLEL2,
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
  OPTION_RESUME,
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
  maxint_t x;
  int64_t a;
  int option;
  int threads;
  std::string resumeFile;
  std::string backupFile;
  bool time;

  CmdOptions() :
    x(-1),
    a(-1),
    option(OPTION_DELEGLISE_RIVAT),
    threads(get_num_threads()),
    time(false)
  { }

  bool is_resume() const
  {
    return !resumeFile.empty();
  }
};

CmdOptions parseOptions(int, char**);

} // namespace

#endif
