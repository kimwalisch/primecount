#include <int128_t.hpp>
#include <Vector.hpp>
#include <BitSieve240.hpp>
#include <calculator.hpp>
#include <primesieve.hpp>
#include <imath.hpp>
#include <iostream>

using namespace primecount;

class Sieve128bit : public BitSieve240 {
public:
  static void sieve128(uint128_t start, uint128_t stop)
  {
    uint128_t old_start = start;
    if (start % 240)
      start -= start % 240;

    uint128_t dist = (stop - start) + 1;
    std::size_t size = (std::size_t) ceil_div(dist, 240);
    uint64_t sqrt_stop = (uint64_t) isqrt(stop);
    Vector<uint64_t> sieve(size);

    std::fill(sieve.begin(), sieve.end(), ~0ull);
    sieve.front() &= unset_smaller_[old_start % 240];
    sieve.back() &= unset_larger_[stop % 240];
    primesieve::iterator iter(7, sqrt_stop);
    uint64_t prime;

    while ((prime = iter.next_prime()) <= sqrt_stop)
    {
      uint128_t q = (start / prime) + 1;
      uint128_t n = prime * q;
      n += prime * (n % 2 == 0);
      ASSERT(n % 2 == 0);

      uint64_t i = (uint64_t) (n - start);
      uint64_t limit = (uint64_t) (stop - start);

      for (; i <= limit; i += prime * 2)
        sieve[i / 240] &= unset_bit_[i % 240];
    }

    for (std::size_t i = 0; i < sieve.size(); i++)
    {
      uint64_t bits = sieve[i];
      for (; bits; bits &= bits - 1)
      {
        auto bit_index = __builtin_ctzll(bits);
        uint64_t bit_value = bit_values_[bit_index];
        uint128_t prime = start + i * 240 + bit_value;
        std::cout << prime << std::endl;
      }
    }
  }
};

int main(int argc, char** argv)
{
  if (argc < 3)
  {
    std::cerr << "Missing start/stop params!" << std::endl;
    return 0;
  }

  uint128_t start = calculator::eval<uint128_t>(argv[1]);
  uint128_t stop = calculator::eval<uint128_t>(argv[2]);
  Sieve128bit::sieve128(start, stop);

  return 0;
}
