///
/// @file  LoadBalancer.hpp
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef LOADBALANCER_HPP
#define LOADBALANCER_HPP

#include <primecount-internal.hpp>
#include <int128_t.hpp>
#include <json.hpp>
#include <Status.hpp>

#include <stdint.h>

namespace primecount {

struct ThreadSettings
{
  int64_t low = 0;
  int64_t segments = 0;
  int64_t segment_size = 0;
  maxint_t sum = 0;
  double init_secs = 0;
  double secs = 0;
  int thread_id = -1;

  void start_time() { secs = get_time(); }
  void stop_time() { secs = get_time() - secs; }
  void init_finished() { init_secs = get_time() - secs; }
};

class LoadBalancer
{
public:
  LoadBalancer(maxint_t x, int64_t y, int64_t z, int64_t k, int64_t sieve_limit, maxint_t sum_approx);
  bool get_work(ThreadSettings& thread);
  void update_result(ThreadSettings& thread);
  void finish_backup();
  bool resume(ThreadSettings& thread) const;
  bool resume(maxint_t& sum, double& time) const;
  int get_resume_threads() const;
  double get_wtime() const;
  maxint_t get_sum() const;

private:
  void backup(ThreadSettings& thread);
  void backup(int thread_id);
  void update(ThreadSettings& thread);
  void update_segments(ThreadSettings& thread);
  double remaining_secs() const;

  maxint_t x_;
  int64_t y_;
  int64_t z_;
  int64_t k_;
  int64_t x_star_;
  int64_t low_;
  int64_t max_low_;
  int64_t sieve_limit_;
  int64_t segments_;
  int64_t segment_size_;
  int64_t max_size_;
  maxint_t sum_;
  maxint_t sum_approx_;
  double time_;
  double backup_time_;
  Status status_;
  nlohmann::json json_;
  nlohmann::json copy_;
};

} // namespace

#endif
