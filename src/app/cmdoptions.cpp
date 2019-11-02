///
/// @file   cmdoptions.cpp
/// @brief  Parse command-line options for the primecount console
///         (terminal) application.
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "cmdoptions.hpp"

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <print.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <cstddef>
#include <map>
#include <string>
#include <type_traits>
#include <vector>
#include <utility>

using namespace std;

namespace primecount {

void help();
void version();
void test();

/// Some command-line options require an additional parameter.
/// Examples: --threads THREADS, -a ALPHA, ...
enum IsParam
{
  NO_PARAM,
  REQUIRED_PARAM,
  OPTIONAL_PARAM
};

/// Command-line options
map<string, std::pair<OptionID, IsParam>> optionMap =
{
  { "-a", make_pair(OPTION_ALPHA, REQUIRED_PARAM) },
  { "--alpha", make_pair(OPTION_ALPHA, REQUIRED_PARAM) },
  { "--alpha-y", make_pair(OPTION_ALPHA_Y, REQUIRED_PARAM) },
  { "--alpha-z", make_pair(OPTION_ALPHA_Z, REQUIRED_PARAM) },
  { "-d", make_pair(OPTION_DELEGLISE_RIVAT, NO_PARAM) },
  { "--deleglise-rivat", make_pair(OPTION_DELEGLISE_RIVAT, NO_PARAM) },
  { "--deleglise-rivat-64", make_pair(OPTION_DELEGLISE_RIVAT_64, NO_PARAM) },
  { "--deleglise-rivat-128", make_pair(OPTION_DELEGLISE_RIVAT_128, NO_PARAM) },
  { "-g", make_pair(OPTION_GOURDON, NO_PARAM) },
  { "--gourdon", make_pair(OPTION_GOURDON, NO_PARAM) },
  { "--gourdon-64", make_pair(OPTION_GOURDON_64, NO_PARAM) },
  { "--gourdon-128", make_pair(OPTION_GOURDON_128, NO_PARAM) },
  { "-h", make_pair(OPTION_HELP, NO_PARAM) },
  { "--help", make_pair(OPTION_HELP, NO_PARAM) },
  { "-l", make_pair(OPTION_LEGENDRE, NO_PARAM) },
  { "--legendre", make_pair(OPTION_LEGENDRE, NO_PARAM) },
  { "--lehmer", make_pair(OPTION_LEHMER, NO_PARAM) },
  { "--lmo", make_pair(OPTION_LMO, NO_PARAM) },
  { "--lmo1", make_pair(OPTION_LMO1, NO_PARAM) },
  { "--lmo2", make_pair(OPTION_LMO2, NO_PARAM) },
  { "--lmo3", make_pair(OPTION_LMO3, NO_PARAM) },
  { "--lmo4", make_pair(OPTION_LMO4, NO_PARAM) },
  { "--lmo5", make_pair(OPTION_LMO5, NO_PARAM) },
  { "-m", make_pair(OPTION_MEISSEL, NO_PARAM) },
  { "--meissel", make_pair(OPTION_MEISSEL, NO_PARAM) },
  { "-n", make_pair(OPTION_NTHPRIME, NO_PARAM) },
  { "--nth-prime", make_pair(OPTION_NTHPRIME, NO_PARAM) },
  { "--number", make_pair(OPTION_NUMBER, REQUIRED_PARAM) },
  { "-p", make_pair(OPTION_PRIMESIEVE, NO_PARAM) },
  { "--primesieve", make_pair(OPTION_PRIMESIEVE, NO_PARAM) },
  { "--Li", make_pair(OPTION_LI, NO_PARAM) },
  { "--Li-inverse", make_pair(OPTION_LIINV, NO_PARAM) },
  { "--Ri", make_pair(OPTION_RI, NO_PARAM) },
  { "--Ri-inverse", make_pair(OPTION_RIINV, NO_PARAM) },
  { "--phi", make_pair(OPTION_PHI, NO_PARAM) },
  { "--P2", make_pair(OPTION_P2, NO_PARAM) },
  { "--S1", make_pair(OPTION_S1, NO_PARAM) },
  { "--S2-easy", make_pair(OPTION_S2_EASY, NO_PARAM) },
  { "--S2-hard", make_pair(OPTION_S2_HARD, NO_PARAM) },
  { "--S2-trivial", make_pair(OPTION_S2_TRIVIAL, NO_PARAM) },
  { "--AC", make_pair(OPTION_AC, NO_PARAM) },
  { "--B", make_pair(OPTION_B, NO_PARAM) },
  { "--D", make_pair(OPTION_D, NO_PARAM) },
  { "--Phi0", make_pair(OPTION_PHI0, NO_PARAM) },
  { "--Sigma", make_pair(OPTION_SIGMA, NO_PARAM) },
  { "-s", make_pair(OPTION_STATUS, OPTIONAL_PARAM) },
  { "--status", make_pair(OPTION_STATUS, OPTIONAL_PARAM) },
  { "--test", make_pair(OPTION_TEST, NO_PARAM) },
  { "--time", make_pair(OPTION_TIME, NO_PARAM) },
  { "-t", make_pair(OPTION_THREADS, REQUIRED_PARAM) },
  { "--threads", make_pair(OPTION_THREADS, REQUIRED_PARAM) },
  { "-v", make_pair(OPTION_VERSION, NO_PARAM) },
  { "--version", make_pair(OPTION_VERSION, NO_PARAM) }
};

/// Command-line option
struct Option
{
  // Example:
  // str = "--threads=32"
  // opt = "--threads"
  // val = "32"
  string str;
  string opt;
  string val;

  template <typename T>
  T to() const
  {
    try {
      if (std::is_floating_point<T>::value)
        return (T) stod(val);
      else
        return (T) to_maxint(val);
    }
    catch (std::exception&) {
      throw primecount_error("invalid option '" + opt + "=" + val + "'");
    }
  }
};

void optionStatus(Option& opt,
                  CmdOptions& opts)
{
  set_print(true);
  opts.time = true;

  if (!opt.val.empty())
    set_status_precision(opt.to<int>());
}

/// e.g. "--thread=4" -> return "--thread"
string getOption(const string& str)
{
  size_t pos = str.find_first_of("=0123456789");

  if (pos == string::npos)
    return str;
  else
    return str.substr(0, pos);
}

/// e.g. "--thread=4" -> return "4"
string getValue(const string& str)
{
  size_t pos = str.find_first_of("0123456789");

  if (pos == string::npos)
    return string();
  else
    return str.substr(pos);
}

/// Parse the next command-line option.
/// e.g. "--threads=32"
/// -> opt.str = "--threads=32"
/// -> opt.opt = "--threads"
/// -> opt.val = "8"
///
Option parseOption(int argc, char* argv[], int& i)
{
  Option opt;
  opt.str = argv[i];

  // Check if the option has the format:
  // --arg or -a (but not --arg=N)
  if (optionMap.count(opt.str))
  {
    opt.opt = opt.str;
    IsParam isParam = optionMap[opt.str].second;

    if (isParam == REQUIRED_PARAM)
    {
      i += 1;

      if (i < argc)
        opt.val = argv[i];
      if (opt.val.empty() || optionMap.count(opt.val))
        throw primecount_error("missing value for option '" + opt.opt + "'");
    }
  }
  else
  {
    // Here the option is either:
    // 1) A number (e.g. the start number)
    // 2) An option of type: --arg=N

    opt.opt = getOption(opt.str);
    opt.val = getValue(opt.str);

    if (opt.opt.empty() && !opt.val.empty())
      opt.opt = "--number";
  }

  if (!optionMap.count(opt.opt))
    throw primecount_error("unrecognized option '" + opt.opt + "'");

  return opt;
}

CmdOptions parseOptions(int argc, char* argv[])
{
  CmdOptions opts;
  vector<maxint_t> numbers;

  for (int i = 1; i < argc; i++)
  {
    Option opt = parseOption(argc, argv, i);
    OptionID optionID = optionMap[opt.opt].first;

    switch (optionID)
    {
      case OPTION_ALPHA:   set_alpha(opt.to<double>()); break;
      case OPTION_ALPHA_Y: set_alpha_y(opt.to<double>()); break;
      case OPTION_ALPHA_Z: set_alpha_z(opt.to<double>()); break;
      case OPTION_NUMBER:  numbers.push_back(opt.to<maxint_t>()); break;
      case OPTION_THREADS: set_num_threads(opt.to<int>()); break;
      case OPTION_HELP:    help(); break;
      case OPTION_STATUS:  optionStatus(opt, opts); break;
      case OPTION_TIME:    opts.time = true; break;
      case OPTION_TEST:    test(); break;
      case OPTION_VERSION: version(); break;
      default:             opts.option = optionID;
    }
  }

  if (opts.option == OPTION_PHI)
  {
    if (numbers.size() >= 2)
      opts.a = numbers[1];
    else
      throw primecount_error("option --phi requires 2 numbers");
  }

  if (numbers.empty())
    throw primecount_error("missing x number");

  opts.x = numbers[0];

  return opts;
}

} // namespace
