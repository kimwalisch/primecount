///
/// @file  backup.cpp
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <backup.hpp>
#include <int128_t.hpp>
#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <WjCryptLib_Md5.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

#if defined(_WIN32) || \
    defined(_WIN64)
  #include <windows.h> // MoveFile
#else
  #include <cstdio> // posix rename file
#endif

using namespace primecount;

namespace {

// Default backup file name
std::string backup_file_ = "primecount.backup";

/// TODO: use C++17 filesystem::rename()
void atomic_backup()
{
#if defined(_WIN32)
  MoveFileEx(std::string(backup_file() + ".new").c_str(), backup_file().c_str(), MOVEFILE_REPLACE_EXISTING);
#else
  // atomically replace old backup file
  std::rename(std::string(backup_file() + ".new").c_str(), backup_file().c_str());
#endif
}

void verify_checksum(nlohmann::json copy)
{
  if (copy.empty())
    return;

  std::string oldMd5 = copy["md5"];
  copy.erase("md5");

  MD5_HASH md5Hash;
  std::string json_str = copy.dump();
  Md5Calculate(json_str.c_str(), (uint32_t) json_str.size(), &md5Hash);

  // convert char buffer into hex std::string
  // https://stackoverflow.com/a/7639754
  std::stringstream md5;
  md5 << std::hex << std::setfill('0');
  for (std::size_t i = 0; i < sizeof(md5Hash.bytes); i++)
    md5 << std::setw(2) << (unsigned) md5Hash.bytes[i];

  if (oldMd5 != md5.str())
    throw primecount_error("MD5 mismatch: corrupted primecount.backup file!");
}

} // namespace

namespace primecount
{

std::string backup_file()
{
  return backup_file_;
}

void set_backup_file(const std::string& filename)
{
  backup_file_ = filename;
}

nlohmann::json load_backup()
{
  std::ifstream ifs(backup_file());
  nlohmann::json j;

  if (ifs.is_open())
    ifs >> j;

  verify_checksum(j);

  return j;
}

void store_backup(nlohmann::json& j)
{
  j.erase("md5");
  MD5_HASH md5Hash;
  std::string json_str = j.dump();
  Md5Calculate(json_str.c_str(), (uint32_t) json_str.size(), &md5Hash);

  // convert char buffer into hex std::string
  // https://stackoverflow.com/a/7639754
  std::stringstream md5;
  md5 << std::hex << std::setfill('0');
  for (std::size_t i = 0; i < sizeof(md5Hash.bytes); i++)
    md5 << std::setw(2) << (unsigned) md5Hash.bytes[i];

  j["md5"] = md5.str();

  {
    // create new backup file
    std::ofstream ofs(backup_file() + ".new");
    ofs << std::setw(4) << j << std::endl;
  }

  atomic_backup();
}

bool is_resume(const nlohmann::json& j,
               const std::string& formula,
               maxint_t x,
               int64_t y)
{
  auto iter = j.find(formula);

  if (iter != j.end())
  {
    auto& jsonFormula = *iter;
    return x == to_maxint(jsonFormula["x"]) &&
           y == jsonFormula["y"]; 
  }

  return false;
}

bool is_resume(const nlohmann::json& j,
               const std::string& formula,
               maxint_t x,
               int64_t y,
               int64_t z,
               int64_t k)
{
  auto iter = j.find(formula);

  if (iter != j.end())
  {
    auto& jsonFormula = *iter;
    return x == to_maxint(jsonFormula["x"]) &&
           y == jsonFormula["y"] &&
           z == jsonFormula["z"] &&
           k == jsonFormula["k"];
  }

  return false;
}

bool is_resume(const nlohmann::json& j,
               const std::string& formula,
               int thread_id,
               maxint_t x,
               int64_t y,
               int64_t z,
               int64_t k)
{
  auto iter = j.find(formula);

  if (iter != j.end())
  {
    auto& jsonFormula = *iter;
    return x == to_maxint(jsonFormula["x"]) &&
           y == jsonFormula["y"] &&
           z == jsonFormula["z"] &&
           k == jsonFormula["k"] &&
           jsonFormula.count("thread" + std::to_string(thread_id)) > 0;
  }

  return false;
}

int calculate_resume_threads(const nlohmann::json& j,
                             const std::string& formula)
{
  auto iter = j.find(formula);
  if (iter == j.end())
    return 0;

  std::string thread = "thread";
  int maxThreadId = -1;
  auto& jsonFormula = *iter;

  for (auto it = jsonFormula.begin(); it != jsonFormula.end(); ++it)
  {
    std::string threadId = it.key();

    // Check if json entry starts with "thread"
    if (!threadId.compare(0, thread.size(), thread))
    {
      // Get the integer part,
      // string format is: "thread<THREAD_ID>"
      threadId.erase(0, thread.size());
      int tid = std::stoi(threadId);
      // std::max<int> prevents collison with max macro from Windows.h
      maxThreadId = std::max<int>(tid, maxThreadId);
    }
  }

  return maxThreadId + 1;
}

} // namespace
