///
/// @file  backup.hpp
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef BACKUP_HPP
#define BACKUP_HPP

#include <string>

namespace primecount
{

std::string backup_file();

void set_backup_file(const std::string& backup_file);

void atomic_backup();

}

#endif
