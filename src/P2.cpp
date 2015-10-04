///
/// @file  P2.cpp
/// @brief 2nd partial sieve function.
///        P2(x, y) counts the numbers <= x that have exactly 2 prime
///        factors each exceeding the a-th prime.
///
/// Copyright (C) 2015 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <aligned_vector.hpp>
#include <BitSieve.hpp>
#include <generate.hpp>
#include <int128.hpp>
#include <min_max.hpp>
#include <pmath.hpp>
#include <Wheel.hpp>

#include <stdint.h>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;

namespace {
namespace P2 {

class ReversePrimeIterator
{
public:
  ReversePrimeIterator(int64_t stop, int64_t start) :
    iter_(stop, start),
    prime_(stop)
  { }
  int64_t previous_prime()
  {
    if (prime_ <= 2)
      return -1;
    prime_ = iter_.previous_prime();
    return prime_;
  }
private:
  primesieve::iterator iter_;
  int64_t prime_;
};

/// Cross-off the multiples inside [low, high[
/// of the primes <= sqrt(high - 1).
///
void cross_off(BitSieve& sieve,
               vector<int32_t> primes,
               Wheel& wheel,
               int64_t c,
               int64_t low,
               int64_t high)
{
  int64_t pi_sqrt_high = pi_bsearch(primes, isqrt(high - 1));

  for (int64_t i = c + 1; i <= pi_sqrt_high; i++)
  {
    int64_t prime = primes[i];
    int64_t m = wheel[i].next_multiple;
    int64_t wheel_index = wheel[i].wheel_index;

    // cross-off the multiples of prime inside [low, high[
    for (; m < high; m += prime * Wheel::next_multiple_factor(&wheel_index))
      sieve.unset(m - low);

    wheel[i].set(m, wheel_index);
  }
}

/// Calculate the segments per thread.
/// The idea is to gradually increase the segments per thread (based
/// on elapsed time) in order to keep all CPU cores busy.
///
int64_t balanceLoad(int64_t segments_per_thread, double start_time)
{
  double seconds = get_wtime() - start_time;

  if (seconds < 30)
    segments_per_thread *= 2;
  else if (segments_per_thread >= 4)
    segments_per_thread -= segments_per_thread / 4;

  return segments_per_thread;
}

template <typename T>
T P2_thread(T x,
            int64_t y,
            int64_t segment_size,
            int64_t segments_per_thread,
            int64_t thread_num,
            int64_t low,
            int64_t limit,
            int64_t& pix,
            int64_t& pix_count,
            vector<int32_t>& primes)
{
  pix = 0;
  pix_count = 0;
  low += thread_num * segments_per_thread * segment_size;
  limit = min(low + segments_per_thread * segment_size, limit);
  int64_t size = pi_bsearch(primes, isqrt(limit)) + 1;
  int64_t start = (int64_t) max(x / limit + 1, y);
  int64_t stop  = (int64_t) min(x / low, isqrt(x));
  T P2_thread = 0;

  // P2_thread = \sum_{i=pi[start]}^{pi[stop]} pi(x / primes[i]) - pi(low - 1)
  // We use a reverse prime iterator to calculate P2_thread
  ReversePrimeIterator prime_iter(stop + 1, start);
  int64_t prime = prime_iter.previous_prime();
  int64_t xp = (int64_t) (x / prime);

  bool sieve_primes = true;
  Wheel wheel(primes, size, low, sieve_primes);
  BitSieve sieve(segment_size);

  // segmented sieve of Eratosthenes
  for (; low < limit; low += segment_size)
  {
    // current segment = interval [low, high[
    int64_t high = min(low + segment_size, limit);
    int64_t c = 6;
    int64_t j = 0;

    // pre-sieve the multiples of the first c primes
    sieve.pre_sieve(c, low, sieve_primes);

    // cross-off the multiples of the primes <= sqrt(high - 1)
    cross_off(sieve, primes, wheel, c, low, high);

    while (prime >= start && 
           xp < high)
    {
      pix += sieve.count(j, xp - low);
      j = xp - low + 1;
      pix_count++;
      P2_thread += pix;
      prime = prime_iter.previous_prime();
      xp = (int64_t) (x / prime);
    }

    pix += sieve.count(j, (high - 1) - low);
  }

  return P2_thread;
}

maxint_t get_next_line(ifstream& infile)
{
  string line;
  getline(infile, line);
  size_t pos = line.find(" = ") + 3;
  return to_maxint(line.substr(pos, line.size() - pos));
}

double get_next_double(ifstream& infile)
{
  string line;
  getline(infile, line);
  size_t pos = line.find(" = ") + 3;
  stringstream ss;
  double d = 0;
  ss << line.substr(pos, line.size() - pos);
  ss >> d;
  return d;
}

template <typename T>
void save_file(T x,
               int64_t y,
               int64_t low,
               int64_t limit,
               int64_t segments_per_thread,
               T pix,
               T p2,
               double time)
{
  ofstream outfile("P2.txt");

  if (outfile.is_open())
  {
    outfile << "x = " << x << endl;
    outfile << "y = " << y << endl;
    outfile << "low = " << low << endl;
    outfile << "limit = " << limit << endl;
    outfile << "segments_per_thread = " << segments_per_thread << endl;
    outfile << "pix = " << pix << endl;
    outfile << "p2 = " << p2 << endl;
    outfile << "Seconds = " << fixed << setprecision(3) << (get_wtime() - time) << endl;
    outfile.close();
  }
}

template <typename T>
void read_file(T x,
               int64_t y,
               int64_t* low,
               int64_t limit,
               int64_t* segments_per_thread,
               T* pix,
               T* p2,
               double* time)
{
  ifstream infile("P2.txt");

  if (infile.is_open())
  {
    try
    {
      T x2 = get_next_line(infile);
      int64_t y2 = (int64_t) get_next_line(infile);
      int64_t low2 = (int64_t) get_next_line(infile);
      int64_t limit2 = (int64_t) get_next_line(infile);
      int64_t segments_per_thread2 = (int64_t) get_next_line(infile);
      T pix2 = get_next_line(infile);
      T p22 = get_next_line(infile);
      double seconds = get_next_double(infile);
      infile.close();

      // only resume if P2.txt matches the
      // command-line values x and alpha
      if (x == x2 &&
          y == y2 &&
          low2 > *low &&
          low2 < limit2 &&
          limit == limit2)
      {
        *low = low2;
        *segments_per_thread = segments_per_thread2;
        *pix = pix2;
        *p2 = p22;
        *time -= seconds;

        if (print_status())
        {
          cout << "=== Resuming from P2.txt ===" << endl;
          cout << "low = " << *low << endl;
          cout << "segments_per_thread = " << *segments_per_thread << endl;
          cout << "pix = " << *pix << endl;
          cout << "p2 = " << *p2 << endl;
          cout << "Seconds = " << seconds << endl;
          cout << endl;
        }
      }
    }
    catch (std::exception&)
    {
      throw primecount_error("failed to read P2.txt");
    }
  }
}

/// P2(x, y) counts the numbers <= x that have exactly 2 prime
/// factors each exceeding the a-th prime, a = pi(y).
/// Space complexity: O((x / y)^(1/2)).
///
template <typename T>
T P2(T x, int64_t y, int threads, double& time)
{
#if __cplusplus >= 201103L
  static_assert(prt::is_signed<T>::value,
                "P2(T x, ...): T must be signed integer type");
#endif

  T a = pi_legendre(y, 1);
  T b = pi_legendre((int64_t) isqrt(x), 1);

  if (x < 4 || a >= b)
    return 0;

  int64_t low = 2;
  int64_t limit = (int64_t)(x / max(y, 1));
  int64_t segment_size = max(isqrt(limit), 1 << 12);
  int64_t segments_per_thread = 64;
  threads = validate_threads(threads, limit);

  vector<int32_t> primes = generate_primes(isqrt(limit));
  aligned_vector<int64_t> pix(threads);
  aligned_vector<int64_t> pix_counts(threads);
  double backup_time = get_wtime();

  // \sum_{i=a+1}^{b} pi(x / primes[i]) - (i - 1)
  T p2 = 0;
  T pix_total = 0;

  // \sum_{i=a+1}^{b} -(i - 1)
  p2 = (a - 2) * (a + 1) / 2 - (b - 2) * (b + 1) / 2;

  read_file(x, y, &low, limit, &segments_per_thread, &pix_total, &p2, &time);

  // \sum_{i=a+1}^{b} pi(x / primes[i])
  while (low < limit)
  {
    int64_t segments = ceil_div(limit - low, segment_size);
    threads = in_between(1, threads, segments);
    segments_per_thread = in_between(1, segments_per_thread, ceil_div(segments, threads));
    double seconds = get_wtime();

    #pragma omp parallel for \
        num_threads(threads) reduction(+: p2)
    for (int i = 0; i < threads; i++)
      p2 += P2_thread(x, y, segment_size, segments_per_thread, i,
         low, limit, pix[i], pix_counts[i], primes);

    low += segments_per_thread * threads * segment_size;
    segments_per_thread = balanceLoad(segments_per_thread, seconds);

    // Add missing sum contributions in order
    for (int i = 0; i < threads; i++)
    {
      p2 += pix_total * pix_counts[i];
      pix_total += pix[i];
    }

    if (get_wtime() - backup_time > 3600)
    {
      save_file(x, y, low, limit, segments_per_thread, pix_total, p2, time);
      backup_time = get_wtime();
    }

    if (print_status())
    {
      int precision = get_status_precision(x);
      double percent = get_percent((double) low, (double) limit);
      cout << "\rStatus: " << fixed << setprecision(precision) << percent << '%' << flush;
    }
  }

  save_file(x, y, limit, limit, segments_per_thread, pix_total, p2, time);

  return p2;
}

} // namespace P2
} // namespace

namespace primecount {

int64_t P2(int64_t x, int64_t y, int threads)
{
  print("");
  print("=== P2(x, y) ===");
  print("Computation of the 2nd partial sieve function");
  print(x, y, threads);

  double time = get_wtime();
  int64_t p2 = P2::P2(x, y, threads, time);

  print("P2", p2, time);
  return p2;
}

#ifdef HAVE_INT128_T

int128_t P2(int128_t x, int64_t y, int threads)
{
  print("");
  print("=== P2(x, y) ===");
  print("Computation of the 2nd partial sieve function");
  print(x, y, threads);

  double time = get_wtime();
  int128_t p2 = P2::P2(x, y, threads, time);

  print("P2", p2, time);
  return p2;
}

#endif

} // namespace primecount
