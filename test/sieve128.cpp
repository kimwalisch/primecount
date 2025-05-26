#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <BitSieve240.hpp>
#include <ctz.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <popcnt.hpp>
#include <Vector.hpp>

#include <iostream>

namespace {

using namespace primecount;

class Sieve_128bit : public BitSieve240
{
public:
  uint64_t count_primes()
  {
    primes_ = 0;
    for (uint64_t bits : sieve_)
      primes_ += popcnt64(bits);

    return primes_;
  }

  /// Sieve interval [low, high]
  template <typename T>
  void sieve(T low, T high)
  {
    T old_low = low;
    if (low % 240)
      low -= low % 240;

    T dist = (high - low) + 1;
    uint64_t size = (uint64_t) ceil_div(dist, 240);
    uint64_t sqrt_high = (uint64_t) isqrt(high);
    uint64_t prime;

    low_ = low;
    sieve_.resize(size);

    std::fill(sieve_.begin(), sieve_.end(), ~0ull);
    sieve_.front() &= unset_smaller_[old_low % 240];
    sieve_.back() &= unset_larger_[high % 240];
    primesieve::iterator iter(7, sqrt_high);

    while ((prime = iter.next_prime()) <= sqrt_high)
    {
      // Calculate first multiple > low
      T q = (low / prime) + 1;
      T n = prime * q;
      n += prime * (n % 2 == 0);
      ASSERT(n % 2 != 0);

      uint64_t i = (uint64_t) (n - low);
      uint64_t limit = (uint64_t) (high - low);

      // Cross-off multiples
      for (; i <= limit; i += prime * 2)
        sieve_[i / 240] &= unset_bit_[i % 240];
    }
  }

  maxuint_t find_nth_prime(uint64_t n) const
  {
    ASSERT(n <= primes_);
    uint64_t count = 0;

    for (std::size_t i = 0; i < sieve_.size(); i++)
    {
      uint64_t bits = sieve_[i];
      uint64_t count_bits = popcnt64(bits);

      if (count + count_bits < n)
        count += count_bits;
      else
      {
        for (; bits; bits &= bits - 1)
        {
          if (++count == n)
          {
            uint64_t bit_index = ctz64(bits);
            uint64_t bit_value = bit_values_[bit_index];
            maxuint_t prime = low_ + i * 240 + bit_value;
            return prime;
          }
        }
      }
    }

    throw primecount_error("Failed to find nth prime!");
  }

private:
  maxuint_t low_ = 0;
  uint64_t primes_ = 0;
  Vector<uint64_t> sieve_;
};

} // namespace

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::cerr << "Missing start param!" << std::endl;
    return 0;
  }

  maxuint_t start = to_maxint(argv[1]);
  maxuint_t prime_approx = start + isqrt(start);
  uint64_t n = (uint64_t) ((prime_approx - start) / ilog(start));
  uint64_t segment_size = (uint64_t) iroot<3>(start) * 30;
  uint64_t count_segments = 0;
  uint64_t count = 0;

  for (maxuint_t low = start; true; low += segment_size)
  {
    count_segments++;
    maxuint_t high = low + segment_size - 1;
    Sieve_128bit sieve;

    if (low <= pstd::numeric_limits<uint64_t>::max() &&
        high <= pstd::numeric_limits<uint64_t>::max())
      sieve.sieve((uint64_t) low, (uint64_t) high);
    else
      sieve.sieve(low, high);

    uint64_t primes = sieve.count_primes();

    if (count + primes < n)
      count += primes;
    else
    {
      maxuint_t nth_prime = sieve.find_nth_prime(n - count);
      
      std::cout << "start: " << start << std::endl;
      std::cout << "segments: " << count_segments << std::endl;
      std::cout << "n: " << n << std::endl;
      std::cout << "nth_prime: " << nth_prime << std::endl;
      return 0;
    }
  }

  return 0;
}
