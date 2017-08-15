///
/// @file  backup.cpp
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <backup.hpp>
#include <string>

#if defined(_WIN32)
  #include <windows.h> // MoveFile
#else
  #include <cstdio> // posix rename file
#endif

using namespace std;

namespace {

string backup_file_ = "primecount.backup";

}

namespace primecount
{

string backup_file()
{
  return backup_file_;
}

void set_backup_file(const string& filename)
{
  backup_file_ = filename;
}

void atomic_backup()
{
  // TODO: use C++17 filesystem::rename()

#if defined(_WIN32)
  MoveFileEx(string(backup_file() + ".new").c_str(), backup_file().c_str(), MOVEFILE_REPLACE_EXISTING);
#else
  // atomically replace old backup file
  rename(string(backup_file() + ".new").c_str(), backup_file().c_str());
#endif
}

} // namespace
