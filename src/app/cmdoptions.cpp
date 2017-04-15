///
/// @file   cmdoptions.cpp
/// @brief  Parse command-line options for the primecount console
///         (terminal) application.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "cmdoptions.hpp"
#include <primecount-internal.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <vector>
#include <string>
#include <map>
#include <exception>
#include <cstdlib>
#include <cstddef>

using std::string;

namespace primecount {

void help();
void version();
bool test();

/// Command-line options
std::map<string, OptionID> optionMap =
{
  { "-a", OPTION_ALPHA },
  { "--alpha", OPTION_ALPHA },
  { "-d", OPTION_DELEGLISE_RIVAT },
  { "--deleglise_rivat", OPTION_DELEGLISE_RIVAT },
  { "--deleglise_rivat1", OPTION_DELEGLISE_RIVAT1 },
  { "--deleglise_rivat2", OPTION_DELEGLISE_RIVAT2 },
  { "--deleglise_rivat_parallel1", OPTION_DELEGLISE_RIVAT_PARALLEL1 },
  { "--deleglise_rivat_parallel2", OPTION_DELEGLISE_RIVAT_PARALLEL2 },
  { "--deleglise_rivat_parallel3", OPTION_DELEGLISE_RIVAT_PARALLEL3 },
  { "-h", OPTION_HELP },
  { "--help", OPTION_HELP },
  { "--legendre", OPTION_LEGENDRE },
  { "--lehmer", OPTION_LEHMER },
  { "-l", OPTION_LMO },
  { "--lmo", OPTION_LMO },
  { "--lmo1", OPTION_LMO1 },
  { "--lmo2", OPTION_LMO2 },
  { "--lmo3", OPTION_LMO3 },
  { "--lmo4", OPTION_LMO4 },
  { "--lmo5", OPTION_LMO5 },
  { "--lmo_parallel1", OPTION_LMO_PARALLEL1 },
  { "--lmo_parallel2", OPTION_LMO_PARALLEL2 },
  { "--lmo_parallel3", OPTION_LMO_PARALLEL3 },
  { "--Li", OPTION_LI },
  { "--Li_inverse", OPTION_LIINV },
  { "-m", OPTION_MEISSEL },
  { "--meissel", OPTION_MEISSEL },
  { "-n", OPTION_NTHPRIME },
  { "--nthprime", OPTION_NTHPRIME },
  { "--number", OPTION_NUMBER },
  { "--P2", OPTION_P2 },
  { "--pi", OPTION_PI },
  { "-p", OPTION_PRIMESIEVE },
  { "--primesieve", OPTION_PRIMESIEVE },
  { "--S1", OPTION_S1 },
  { "--S2_easy", OPTION_S2_EASY },
  { "--S2_hard", OPTION_S2_HARD },
  { "--S2_trivial", OPTION_S2_TRIVIAL },
  { "-s", OPTION_STATUS },
  { "--status", OPTION_STATUS },
  { "--test", OPTION_TEST },
  { "--time", OPTION_TIME },
  { "-t", OPTION_THREADS },
  { "--threads", OPTION_THREADS },
  { "-v", OPTION_VERSION },
  { "--version", OPTION_VERSION }
};

/// e.g. str = "--threads", value = "4"
struct Option
{
  string str;
  string value;
  template <typename T>
  T getValue() const
  {
    return (T) to_maxint(value);
  }
};

void optionStatus(Option& opt,
                  CmdOptions& opts)
{
  set_print_status(true);
  opts.time = true;

  if (!opt.value.empty())
    set_status_precision(opt.getValue<int>());
}

/// e.g. "--threads=8"
/// -> opt.str = "--threads"
/// -> opt.value = "8"
///
Option makeOption(const string& str)
{
  Option opt;
  size_t delimiter = string::npos;
  if (optionMap.count(str) == 0)
    delimiter = str.find_first_of("=0123456789");

  if (delimiter == string::npos)
    opt.str = str;
  else
  {
    opt.str = str.substr(0, delimiter);
    opt.value = str.substr(delimiter + (str.at(delimiter) == '=' ? 1 : 0));
  }
  if (opt.str.empty() && !opt.value.empty())
    opt.str = "--number";
  if (!optionMap.count(opt.str))
    opt.str = "--help";

  return opt;
}

CmdOptions parseOptions(int argc, char* argv[])
{
  CmdOptions opts;
  std::vector<maxint_t> numbers;

  try
  {
    for (int i = 1; i < argc; i++)
    {
      Option opt = makeOption(argv[i]);

      switch (optionMap[opt.str])
      {
        case OPTION_ALPHA:   set_alpha(std::stod(opt.value)); break;
        case OPTION_NUMBER:  numbers.push_back(opt.getValue<maxint_t>()); break;
        case OPTION_THREADS: opts.threads = opt.getValue<int>(); break;
        case OPTION_HELP:    help(); break;
        case OPTION_STATUS:  optionStatus(opt, opts) ;break;
        case OPTION_TIME:    opts.time = true; break;
        case OPTION_TEST:    if (test()) exit(0); exit(1);
        case OPTION_VERSION: version(); break;
        default:             opts.option = optionMap[opt.str];
      }
    }
  }
  catch (std::exception&)
  {
    help();
  }

  if (numbers.size() == 1)
    opts.x = numbers[0];
  else
    help();

  return opts;
}

} // namespace
