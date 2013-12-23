///
/// @file  phi.cpp
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "pi_bsearch.h"
#include "imath.h"

#include <primesieve.hpp>
#include <stdint.h>
#include <vector>
#include <limits>
#include <cstddef>
#include <cassert>

#ifdef _OPENMP
  #include <omp.h>
  #include "to_omp_threads.h"
#endif

/// Keep the cache size below CACHE_BYTES_LIMIT per thread
#define CACHE_BYTES_LIMIT (16 << 20)

/// Cache phi(x, a) results if a <= CACHE_A_LIMIT
#define CACHE_A_LIMIT 500

/// Avoid slow 64-bit division if possible
#define FAST_DIV(x, y) ((x <= std::numeric_limits<uint32_t>::max()) \
    ? static_cast<uint32_t>(x) / (y) : (x) / (y))

namespace primecount {

/// This class calculates phi(x, a) in constant time for small values
/// of a < 7 using lookup tables. Let pp = prime_products_[a]:
/// phi(x, a) = (x / pp) * φ(pp) + phi(x % pp, a).
/// 
class PhiSmall {
public:
  PhiSmall()
  {
    phi_cache_[0].push_back(0);

    // Initialize the phi_cache_ lookup tables
    for (int a = 1; is_cached(a); a++)
      for (int x = 0; x < prime_products_[a]; x++)
        phi_cache_[a].push_back(static_cast<int16_t>(phi(x, a - 1) - phi(x / primes_[a], a - 1)));
  }
  bool is_cached(int64_t a) const
  {
    return a >= 0 && a < 7;
  }
  /// Partial sieve function (a.k.a. Legendre-sum).
  /// phi(x, a) counts the numbers <= x that are not divisible
  /// by any of the first a primes.
  /// @pre is_cached(a) == true.
  ///
  int64_t phi(int64_t x, int64_t a) const
  {
    assert(x >= 0);
    assert(is_cached(a));
    return (x / prime_products_[a]) * totients_[a] + phi_cache_[a][x % prime_products_[a]];
  }
private:
  std::vector<int16_t> phi_cache_[7];
  static const int32_t primes_[7];
  static const int32_t prime_products_[7];
  static const int32_t totients_[7];
};

const int32_t PhiSmall::primes_[7] = { -1, 2, 3, 5, 7, 11, 13 };

/// prime_products_[n] = \prod_{i=1}^{n} primes_[i]
const int32_t PhiSmall::prime_products_[7] = { 1, 2, 6, 30, 210, 2310, 30030 };

/// totients_[n] = \prod_{i=1}^{n} (primes_[i] - 1)
const int32_t PhiSmall::totients_[7] = { 1, 1, 2, 8, 48, 480, 5760 };

/// This class calculates phi(x, a) using the recursive formula:
/// phi(x, a) = phi(x, a - 1) - phi(x / primes_[a], a - 1).
/// This implementation is based on an algorithm from Tomas Oliveira e
/// Silva [1]. I have added a cache to my implementation in which
/// results of phi(x, a) are stored if x < 2^16 and a <= 500.
/// The cache speeds up the calculations by at least 3 orders of
/// magnitude near 10^15.
///
/// [1] Tomas Oliveira e Silva, "Computing pi(x): the combinatorial method",
///     Revista do DETUA, vol. 4, no. 6, pp. 759-768, March 2006
///
class PhiCache {
public:
  PhiCache(const std::vector<int32_t>& primes, const PhiSmall& phiSmall) :
    primes_(primes),
    phiSmall_(phiSmall),
    bytes_(0)
  {
    std::size_t max_size = CACHE_A_LIMIT + 1;
    cache_.resize(std::min(primes.size(), max_size));
  }

  /// phi(x, a) = pi(x) - a + 1
  int64_t phi_bsearch(int64_t x, int64_t a)
  {
    int64_t pix = pi_bsearch(primes_.begin(), primes_.end(), x);
    return pix - a + 1;
  }

  template<int64_t SIGN>
  int64_t phi(int64_t x, int64_t a)
  {
    if (phiSmall_.is_cached(a))
      return phiSmall_.phi(x, a) * SIGN;

    int64_t sum = x * SIGN;
    // phi(x, i) = 1 for iters <= i < a
    int64_t iters = pi_bsearch(primes_.begin(), primes_.begin() + a, isqrt(x));
    sum += (a - iters) * -SIGN;

    for (int64_t a2 = 0; a2 < iters; a2++)
    {
      // x2 = x / primes_[a2]
      int64_t x2 = FAST_DIV(x, primes_[a2]);
      int64_t phi_result;

      if (validate(a2, x2) && cache_[a2][x2] != 0)
        phi_result = cache_[a2][x2] * -SIGN;
      else
      {
        if (x2 <= primes_.back() && x2 < isquare(primes_[a2]))
          phi_result = phi_bsearch(x2, a2) * -SIGN;
        else
          phi_result = phi<-SIGN>(x2, a2);

        if (validate(a2, x2))
          cache_[a2][x2] = static_cast<uint16_t>(phi_result * -SIGN);
      }
      sum += phi_result;
    }
    return sum;
  }
private:
  /// Cache for phi(x, a) results
  std::vector<std::vector<uint16_t> > cache_;
  const std::vector<int32_t>& primes_;
  const PhiSmall& phiSmall_;
  int64_t bytes_;

  bool validate(int64_t a2, int64_t x2)
  {
    if (a2 > CACHE_A_LIMIT || x2 > std::numeric_limits<uint16_t>::max())
      return false;
    // resize and initialize cache if necessary
    if (x2 >= static_cast<int64_t>(cache_[a2].size()))
    {
      if (bytes_ > CACHE_BYTES_LIMIT)
        return false;
      bytes_ += (x2 + 1 - cache_[a2].size()) * 2;
      cache_[a2].resize(x2 + 1, 0);
    }
    return true;
  }
};

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
///
int64_t phi(int64_t x, int64_t a, int threads)
{
  if (x < 1) return 0;
  if (a > x) return 1;
  if (a < 1) return x;

  /// @warning This is only thread safe since C++11
  static const PhiSmall phiSmall;

  if (phiSmall.is_cached(a))
    return phiSmall.phi(x, a);

  std::vector<int32_t> primes;
  primesieve::generate_n_primes(a, &primes);

  if (primes.back() >= x)
    return 1;

  int iters = pi_bsearch(primes.begin(), primes.end(), isqrt(x));
  PhiCache cache(primes, phiSmall);
  int64_t sum = x - a + iters;

#ifdef _OPENMP
  threads = to_omp_threads(threads);
  #pragma omp parallel for firstprivate(cache) reduction(+: sum) \
      num_threads(threads) schedule(dynamic, 16)
#endif
  for (int64_t a2 = 0; a2 < iters; a2++)
    sum += cache.phi<-1>(x / primes[a2], a2);

  return sum;
}

} // namespace primecount
