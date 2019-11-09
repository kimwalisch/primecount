///
/// @file  backup.hpp
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef BACKUP_HPP
#define BACKUP_HPP

#include <stdint.h>
#include <string>

#include <json.hpp>
#include <int128_t.hpp>

namespace primecount
{

std::string backup_file();
void set_backup_file(const std::string& backup_file);
nlohmann::json load_backup();
void store_backup(nlohmann::json& j);
bool is_resume(const nlohmann::json& j, const std::string& formula, maxint_t x, int64_t y);
bool is_resume(const nlohmann::json& j, const std::string& formula, maxint_t x, int64_t y, int64_t z, int64_t k);
bool is_resume(const nlohmann::json& j, const std::string& formula, int thread_id, maxint_t x, int64_t y, int64_t z, int64_t k);
int calculate_resume_threads(const nlohmann::json& j, const std::string& formula);

} // namespace

#endif
