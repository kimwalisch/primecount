///
/// @file  StatusS2.hpp
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef STATUSS2_HPP
#define STATUSS2_HPP

#include <int128_t.hpp>

namespace primecount {

/// Models the work distribution of D, S2_hard and LMO
/// to estimate their progress: evenly distributed sieving
/// work plus the special leaves' work, which starts at
/// leaves_start and then decays like a Lomax distribution.
/// Each algorithm uses its own fitted constants, see
/// doc/StatusS2-tuning.md.
///
struct WorkModelS2
{
  WorkModelS2(double leaves_start,
              double sieve_work,
              double leaves_scale,
              double leaves_shape,
              double sieve_exponent)
    : leaves_start(leaves_start),
      sieve_work(sieve_work),
      leaves_scale(leaves_scale),
      leaves_shape(leaves_shape),
      sieve_exponent(sieve_exponent)
  { }

  // Where the work of the special leaves becomes significant
  double leaves_start;
  // Sieving work relative to the special leaves' work, alpha-adjusted
  double sieve_work;
  // Lomax scale relative to leaves_start
  double leaves_scale;
  // Lomax shape, the leaves' work decays like n^-(1 + leaves_shape)
  double leaves_shape;
  // sieve_work grows like (sieve_limit / leaves_start)^sieve_exponent
  double sieve_exponent;
};

class StatusS2
{
public:
  StatusS2(maxint_t x, int threads, bool is_print);
  StatusS2(maxint_t x, int64_t sieve_limit, const WorkModelS2& model, int threads, bool is_print);
  double getPercent(int64_t low) const;
  void print_S2_hard(double percent);
  void print_S2_easy(int64_t b, int64_t max_b);

private:
  double work(double low) const;
  double leaves_start_ = 0;
  double leaves_scale_ = 1;
  double leaves_shape_ = 1;
  double sieve_work_ = 0;
  double sieve_limit_ = 1;
  double total_work_ = 1;
  double epsilon_ = 0;
  double percent_ = -1;
  double time_ = 0;
  // Only print status if 0.1 seconds have elapsed
  // since last printing the status.
  double threshold_ = 0.1;
  int precision_ = 0;
};

} // namespace

#endif
