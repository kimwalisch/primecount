#ifndef NEXT_N_PRIMES_VECTOR_H
#define NEXT_N_PRIMES_VECTOR_H

#include <primesieve/soe/PrimeSieve.h>
#include <primesieve/soe/PrimeSieveCallback.h>
#include <stdint.h>
#include <vector>
#include <exception>

namespace primecount {

template <typename T>
class Next_N_Primes_Vector : public std::vector<T>, 
                             public PrimeSieveCallback<uint64_t> {
public:
  void generatePrimes(uint64_t start, uint64_t n)
  {
    n_ = n;
    try {
      while (n_ > 0)
      {
        uint64_t stop = start + n_ * 50;
        PrimeSieve ps;
        ps.generatePrimes(start, stop, this);
        start = stop + 1;
      }
    } catch (stop_primesieve&) { }
  }

  void callback(uint64_t prime)
  {
    this->push_back( static_cast<T>(prime) );
    if (--n_ == 0)
      throw stop_primesieve();
  }
private:
  uint64_t n_;
  class stop_primesieve : public std::exception { };
};

} // namespace primecount

#endif
