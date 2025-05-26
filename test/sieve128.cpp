#include <primecount.hpp>
#include <primecount-config.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <BitSieve240.hpp>
#include <ctz.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <popcnt.hpp>
#include <Vector.hpp>

#include <algorithm>
#include <iostream>
#include <omp.h>

namespace {

using namespace primecount;

class Sieve_128bit : public BitSieve240
{
public:
  uint64_t get_prime_count() const
  {
    return count_;
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
    count_ = 0;
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
  
    // Count primes (1 bits)
    for (uint64_t bits : sieve_)
      count_ += popcnt64(bits);
  }

  maxuint_t find_nth_prime_forward(uint64_t n) const
  {
    ASSERT(n > 0);
    ASSERT(n <= count_);

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

  maxuint_t find_nth_prime_backward(uint64_t n) const
  {
    ASSERT(n > 0);
    ASSERT(n <= count_);

    uint64_t count = 0;
    int64_t size = (int64_t) sieve_.size();

    for (int64_t i = size - 1; i >= 0; i--)
    {
      uint64_t bits = sieve_[i];
      uint64_t count_bits = popcnt64(bits);

      if (count + count_bits < n)
        count += count_bits;
      else
      {
        Vector<maxuint_t> primes;
        primes.reserve(count_bits);

        for (; bits; bits &= bits - 1)
        {
          uint64_t bit_index = ctz64(bits);
          uint64_t bit_value = bit_values_[bit_index];
          maxuint_t prime = low_ + i * 240 + bit_value;
          primes.push_back(prime);
        }

        uint64_t j = n - count;
        return primes[primes.size() - j];
      }
    }

    throw primecount_error("Failed to find nth prime!");
  }

private:
  maxuint_t low_ = 0;
  uint64_t count_ = 0;
  Vector<uint64_t> sieve_;
};

/// The aligned_vector class aligns each of its
/// elements on a new cache line in order to avoid
/// false sharing (cache trashing) when multiple
/// threads write to adjacent elements.
///
template <typename T>
class aligned_vector
{
  static_assert(sizeof(T) < MAX_CACHE_LINE_SIZE,
                "sizeof(T) must be < MAX_CACHE_LINE_SIZE");

public:
  aligned_vector(std::size_t size) : vect_(size) { }
  std::size_t size() const { return vect_.size(); }
  T& operator[](std::size_t pos) { return vect_[pos].val; }

private:
  struct CacheLine {
    T val;
    // We cannot use alignas(MAX_CACHE_LINE_SIZE) for
    // the CacheLine struct as GCC does not yet support
    // alignas(n) with n > 128. Also alignas(n) for
    // over-aligned data and dynamic memory allocation
    // is only supported since C++17.
    MAYBE_UNUSED char pad[MAX_CACHE_LINE_SIZE - sizeof(T)];
  };

  Vector<CacheLine> vect_;
};

/// Find the nth prime >= start
maxuint_t find_nth_prime_forward(uint64_t n, maxuint_t start)
{
  ASSERT(n > 0);

  maxuint_t nth_prime = 0;
  uint64_t primes = 0;
  uint64_t while_iters = 0;
  uint64_t segment_size = (uint64_t) iroot<3>(start) * 30;
  uint64_t min_segment_size = (uint64_t) 1e7;
  segment_size = std::max(min_segment_size, segment_size);

  uint64_t avg_prime_gap = ilog(start) + 2;
  uint64_t dist_approx = n * avg_prime_gap;

  int threads = get_num_threads();
  threads = ideal_num_threads(dist_approx, threads, segment_size);
  aligned_vector<Sieve_128bit> sieves(threads);
  bool finished = false;

  #pragma omp parallel num_threads(threads)
  while (!finished)
  {
    int thread_id = omp_get_thread_num();
    uint64_t i = while_iters * threads + thread_id;
    maxuint_t low = start + i * segment_size;
    maxuint_t high = low + segment_size - 1;

    if (low <= pstd::numeric_limits<uint64_t>::max() &&
        high <= pstd::numeric_limits<uint64_t>::max())
      sieves[thread_id].sieve((uint64_t) low, (uint64_t) high);
    else
      sieves[thread_id].sieve(low, high);

    // Wait until all threads have finished
    // computing their current segment.
    #pragma omp barrier
    #pragma omp master
    {
      while_iters++;

      for (int j = 0; j < threads; j++)
      {
        if (primes + sieves[j].get_prime_count() < n)
          primes += sieves[j].get_prime_count();
        else
        {
          nth_prime = sieves[j].find_nth_prime_forward(n - primes);
          finished = true;
          break;
        }
      }
    }

    // Other threads must wait until master
    // thread finishes single-threaded section.
    #pragma omp barrier
  }

  return nth_prime;
}

/// Find the nth prime <= start
maxuint_t find_nth_prime_backward(uint64_t n, maxuint_t start)
{
  ASSERT(n > 0);

  maxuint_t nth_prime = 0;
  uint64_t primes = 0;
  uint64_t while_iters = 0;
  uint64_t segment_size = (uint64_t) iroot<3>(start) * 30;
  uint64_t min_segment_size = (uint64_t) 1e7;
  segment_size = std::max(min_segment_size, segment_size);

  uint64_t avg_prime_gap = ilog(start) + 2;
  uint64_t dist_approx = n * avg_prime_gap;
  if (dist_approx > start)
    dist_approx = (uint64_t) start;

  int threads = get_num_threads();
  threads = ideal_num_threads(dist_approx, threads, segment_size);
  aligned_vector<Sieve_128bit> sieves(threads);
  bool finished = false;

  #pragma omp parallel num_threads(threads)
  while (!finished)
  {
    int thread_id = omp_get_thread_num();
    uint64_t i = while_iters * threads + thread_id;

    if (start > i * segment_size)
    {
      maxuint_t high = start - i * segment_size;
      maxuint_t low = 0;

      if (high >= segment_size)
        low = (high - segment_size) + 1;

      if (low <= pstd::numeric_limits<uint64_t>::max() &&
          high <= pstd::numeric_limits<uint64_t>::max())
        sieves[thread_id].sieve((uint64_t) low, (uint64_t) high);
      else
        sieves[thread_id].sieve(low, high);
    }

    // Wait until all threads have finished
    // computing their current segment.
    #pragma omp barrier
    #pragma omp master
    {
      while_iters++;

      for (int j = 0; j < threads; j++)
      {
        if (primes + sieves[j].get_prime_count() < n)
          primes += sieves[j].get_prime_count();
        else
        {
          nth_prime = sieves[j].find_nth_prime_backward(n - primes);
          finished = true;
          break;
        }
      }
    }

    // Other threads must wait until master
    // thread finishes single-threaded section.
    #pragma omp barrier
  }

  return nth_prime;
}

} // namespace

int main(int argc, char** argv)
{
  if (argc < 3)
  {
    std::cerr << "Missing start/n params!" << std::endl;
    return 0;
  }

  maxuint_t start = to_maxint(argv[1]);
  int64_t n = (int64_t) to_maxint(argv[2]);
  maxuint_t nth_prime;

  std::cout << "n: " << n << std::endl;
  std::cout << "start: " << start << std::endl;

  if (n > 0)
    nth_prime = find_nth_prime_forward(n, start);
  else
    nth_prime = find_nth_prime_backward(-n, start);

  std::cout << "nth_prime: " << nth_prime << std::endl;

  return 0;
}
