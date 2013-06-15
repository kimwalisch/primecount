#ifndef PRIMESIEVEVECTOR_H
#define PRIMESIEVEVECTOR_H

#include <primesieve/soe/PrimeSieve.h>
#include <primesieve/soe/PrimeSieveCallback.h>
#include <stdint.h>
#include <vector>
#include <exception>

namespace primecount {

  template <typename T>
  class PrimeSieveVector : public std::vector<T>,
  public PrimeSieveCallback<uint64_t> {
  public:
    void generatePrimes(uint64_t start, uint64_t stop)
    {
      if (stop >= start)
      {
        // huge number > number of primes
        n_ = stop - start + 2;
        PrimeSieve ps;
        ps.generatePrimes(start, stop, this);
      }
    }
    void generate_N_Primes(uint64_t start, uint64_t n)
    {
      n_ = n;
      this->reserve(this->size() + n_);
      try {
        while (n_ > 0)
        {
          uint64_t stop = start + n_ * 50 + 10000;
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
