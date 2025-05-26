#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <BitSieve240.hpp>
#include <ctz.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <Vector.hpp>

#include <iostream>

namespace {

using namespace primecount;

class Sieve128bit : public BitSieve240
{
public:
  template <typename T>
  static void sieve128(T start, T stop)
  {
    T old_start = start;
    if (start % 240)
      start -= start % 240;

    T dist = (stop - start) + 1;
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
      T q = (start / prime) + 1;
      T n = prime * q;
      n += prime * (n % 2 == 0);
      ASSERT(n % 2 != 0);

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
        uint64_t bit_index = ctz64(bits);
        uint64_t bit_value = bit_values_[bit_index];
        T prime = start + i * 240 + bit_value;
        std::cout << prime << std::endl;
      }
    }
  }
};

} // namespace

int main(int argc, char** argv)
{
  if (argc < 3)
  {
    std::cerr << "Missing start/stop params!" << std::endl;
    return 0;
  }

  maxuint_t start = to_maxint(argv[1]);
  maxuint_t stop = to_maxint(argv[2]);

  if (start <= pstd::numeric_limits<uint64_t>::max() &&
      stop <= pstd::numeric_limits<uint64_t>::max())
    Sieve128bit::sieve128((uint64_t) start, (uint64_t) stop);
  else
    Sieve128bit::sieve128(start, stop);

  return 0;
}
