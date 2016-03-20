///
/// @file   LibdividePrimes.hpp
/// @brief  The LibdividePrimes class calculates (x / primes[i])
///         using libdivide if possible. libdivide allows to
///         replace expensive integer divides with comparatively
///         cheap multiplication and bitshifts.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef LIBDIVIDEPRIMES_HPP
#define LIBDIVIDEPRIMES_HPP

#include <libdivide.h>
#include <generate.hpp>

#include <limits>
#include <stdint.h>
#include <vector>
#include <unistd.h>

namespace primecount {

template <typename T>
struct LibdividePrimes
{
public:
  LibdividePrimes(int64_t max)
  {
    primes_ = generate_primes<T>(max);
    fastdiv_.assign(primes_.begin(), primes_.end());
  }

  /// @return primes[i]
  T operator[](size_t i) const
  {
    return primes_[i];
  }

  /// Compute x / primes[i] using libdivide if possible
  /// otherwise use normal integer division.
  ///
  template <typename X>
  X libdivide(X x, uint64_t i) const
  {
    if (x <= std::numeric_limits<uint64_t>::max())
      return ((uint64_t) x) / fastdiv_[i];

    return x / primes_[i];
  }

private:
  std::vector<T> primes_;
  std::vector<libdivide::divider<uint64_t> > fastdiv_;
};

} // namespace

#endif /* LIBDIVIDEPRIMES_HPP */
