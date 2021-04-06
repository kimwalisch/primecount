///
/// @file  LoadBalancerAC.hpp
/// @brief This load balancer assigns work to the threads in the
///        computation of the A & C formulas (AC.cpp) in
///        Xavier Gourdon's algorithm.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef LOADBALANCERAC_HPP
#define LOADBALANCERAC_HPP

#include <OmpLock.hpp>
#include <stdint.h>

namespace primecount {

class LoadBalancerAC
{
public:
  LoadBalancerAC(int64_t sqrtx, int64_t y, int threads);
  bool get_work(int64_t& low, int64_t& high);

private:
  int64_t low_;
  int64_t sqrtx_;
  int64_t x14_;
  int64_t y_;
  int threads_;
  OmpLock lock_;
};

} // namespace

#endif
