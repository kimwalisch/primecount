#ifndef PRIMECOUNTVECTOR_H
#define PRIMECOUNTVECTOR_H

#include "upper_bound.h"

#include <primesieve/soe/PrimeSieveCallback.h>
#include <stdint.h>
#include <vector>
#include <algorithm>
#include <iterator>

namespace primecount {

template <typename T>
class PrimeCountVector : public std::vector<T>,
                         public PrimeSieveCallback<uint64_t> {
public:
  void callback(uint64_t prime)
  {
    this->push_back( static_cast<T>(prime) );
  }
  T pi(T n) const
  {
    return static_cast<T>(
        std::distance(this->begin(),
                      std::upper_bound(this->begin(), this->end(), n))
    );
  }
};
  
} // end namespace

#endif
