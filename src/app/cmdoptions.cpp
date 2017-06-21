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

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <backup.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cstddef>

using namespace std;

namespace primecount {

void help();
void version();
void test();

/// Command-line options
map<string, OptionID> optionMap =
{
  { "-a", OPTION_ALPHA },
  { "--alpha", OPTION_ALPHA },
  { "-b", OPTION_BACKUP },
  { "--backup", OPTION_BACKUP },
  { "-d", OPTION_DELEGLISE_RIVAT },
  { "--deleglise_rivat", OPTION_DELEGLISE_RIVAT },
  { "--deleglise_rivat1", OPTION_DELEGLISE_RIVAT1 },
  { "--deleglise_rivat2", OPTION_DELEGLISE_RIVAT2 },
  { "--deleglise_rivat_parallel1", OPTION_DELEGLISE_RIVAT_PARALLEL1 },
  { "--deleglise_rivat_parallel2", OPTION_DELEGLISE_RIVAT_PARALLEL2 },
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
  { "--lmo_parallel", OPTION_LMO_PARALLEL },
  { "--Li", OPTION_LI },
  { "--Li_inverse", OPTION_LIINV },
  { "-m", OPTION_MEISSEL },
  { "--meissel", OPTION_MEISSEL },
  { "-n", OPTION_NTHPRIME },
  { "--nthprime", OPTION_NTHPRIME },
  { "--number", OPTION_NUMBER },
  { "--P2", OPTION_P2 },
  { "--phi", OPTION_PHI },
  { "--pi", OPTION_PI },
  { "-p", OPTION_PRIMESIEVE },
  { "--primesieve", OPTION_PRIMESIEVE },
  { "-r", OPTION_RESUME },
  { "--resume", OPTION_RESUME },
  { "--Ri", OPTION_RI },
  { "--Ri_inverse", OPTION_RIINV },
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

/// Command-line option
/// e.g. opt = "--threads", val = "4"
struct Option
{
  string str;
  string opt;
  string val;
  template <typename T>
  T getValue() const
  {
    if (val.empty())
      throw primecount_error("missing value for option " + str);
    return (T) to_maxint(val);
  }
};

void optionStatus(Option& opt,
                  CmdOptions& opts)
{
  set_print(true);
  opts.time = true;

  if (!opt.val.empty())
    set_status_precision(opt.getValue<int>());
}

/// e.g. "--thread=4" -> return "--thread"
string getOption(string str)
{
  size_t pos = str.find_first_of("=0123456789");

  if (pos == string::npos)
    return str;
  else
    return str.substr(0, pos);
}

/// e.g. "--thread=4" -> return "4"
string getValue(string str)
{
  size_t pos = str.find_first_of("0123456789");

  if (pos == string::npos)
    return string();
  else
    return str.substr(pos);
}

/// e.g. "--resume=primecount.backup" -> return "primecount.backup"
string getBackupFile(string str)
{
  size_t pos = str.find_first_of("=");

  if (pos != string::npos)
    return str.substr(pos + 1);

  return backup_file();
}

void optionBackup(Option& opt,
                  CmdOptions& opts)
{
  opts.backupFile = getBackupFile(opt.str);
  set_backup_file(opts.backupFile);
}

void optionResume(Option& opt,
                  CmdOptions& opts)
{
  opts.resumeFile = getBackupFile(opt.str);
  set_backup_file(opts.resumeFile);

  ifstream ifs(backup_file());

  if (!ifs.is_open())
    throw primecount_error("failed to open backup file: " + backup_file());
}

/// e.g. "--threads=8"
/// -> opt.opt = "--threads"
/// -> opt.val = "8"
///
Option makeOption(const string& str)
{
  Option opt;
  opt.str = str;

  if (optionMap.count(str))
    opt.opt = str;
  else
  {
    opt.opt = getOption(str);
    opt.val = getValue(str);
  }

  if (opt.opt.empty() && !opt.val.empty())
    opt.opt = "--number";

  if (!optionMap.count(opt.opt))
    throw primecount_error("unknown option " + str);

  return opt;
}

CmdOptions parseOptions(int argc, char* argv[])
{
  CmdOptions opts;
  vector<maxint_t> numbers;

  for (int i = 1; i < argc; i++)
  {
    Option opt = makeOption(argv[i]);

    switch (optionMap[opt.opt])
    {
      case OPTION_ALPHA:   set_alpha(stod(opt.val)); break;
      case OPTION_BACKUP:  optionBackup(opt, opts); break;
      case OPTION_NUMBER:  numbers.push_back(opt.getValue<maxint_t>()); break;
      case OPTION_THREADS: opts.threads = opt.getValue<int>(); break;
      case OPTION_PHI:     opts.a = opt.getValue<int64_t>(); opts.option = OPTION_PHI; break;
      case OPTION_HELP:    help(); break;
      case OPTION_RESUME:  optionResume(opt, opts); break;
      case OPTION_STATUS:  optionStatus(opt, opts); break;
      case OPTION_TIME:    opts.time = true; break;
      case OPTION_TEST:    test(); break;
      case OPTION_VERSION: version(); break;
      default:             opts.option = optionMap[opt.opt];
    }
  }

  if (!opts.is_resume())
  {
    if (numbers.empty())
      throw primecount_error("missing x number");
    else
      opts.x = numbers[0];
  }

  if (!opts.backupFile.empty() &&
      !opts.resumeFile.empty() &&
      opts.backupFile != opts.resumeFile)
  {
    throw primecount_error("resume and backup file must be the same");
  }

  return opts;
}

} // namespace
