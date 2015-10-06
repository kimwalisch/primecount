///
/// @file  S2_trivial.cpp
/// @brief Calculate the contribution of the trivial special leaves
///        in parallel using OpenMP.
///
/// Copyright (C) 2015 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <generate.hpp>
#include <int128.hpp>

#include <stdint.h>
#include <algorithm>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;

namespace {
namespace S2_trivial {

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
               int64_t c,
               T s2_trivial,
               double time)
{
  ofstream outfile("S2_trivial.txt");

  if (outfile.is_open())
  {
    outfile << "x = " << x << endl;
    outfile << "y = " << y << endl;
    outfile << "c = " << c << endl;
    outfile << "S2_trivial = " << s2_trivial << endl;
    outfile << "Seconds = " << fixed << setprecision(3) << (get_wtime() - time) << endl;
    outfile << "Status = 100%" << endl;
    outfile.close();
  }
}

template <typename T>
bool read_file(T x,
               int64_t y,
               int64_t c,
               T* s2_trivial,
               double* time)
{
  ifstream infile("S2_trivial.txt");

  if (infile.is_open())
  {
    try
    {
      T x2 = get_next_line(infile);
      int64_t y2 = (int64_t) get_next_line(infile);
      int64_t c2 = (int64_t) get_next_line(infile);
      int64_t s2_trivial2 = (int64_t) get_next_line(infile);
      double seconds = get_next_double(infile);

      infile.close();

      if (x == x2 &&
          y == y2 &&
          c == c2)
      {
        *s2_trivial = s2_trivial2;
        *time -= seconds;

        if (print_status())
        {
          if (!print_variables())
            cout << endl;

          cout << "--- Resuming from S2_trivial.txt ---" << endl;
          cout << "x = " << x << endl;
          cout << "y = " << y << endl;
          cout << "c = " << c << endl;
          cout << "S2_trivial = " << *s2_trivial << endl;
          cout << "Seconds = " << seconds << endl;
          cout << "Status = 100%" << endl;
          cout << endl;

          return true;
        }
      }
    }
    catch (std::exception&)
    {
      throw primecount_error("failed to read S2_trivial.txt");
    }
  }

  return false;
}

template <typename T>
T S2_trivial(T x,
             int64_t y,
             int64_t z,
             int64_t c,
             int threads)
{
  print("");
  print("=== S2_trivial(x, y) ===");
  print("Computation of the trivial special leaves");
  print(x, y, c, threads);

  T s2_trivial = 0;
  double time = get_wtime();

  if(!read_file(x, y, c, &s2_trivial, &time))
  {
    int64_t thread_threshold = ipow(10, 7);
    threads = validate_threads(threads, y, thread_threshold);

    PiTable pi(y);
    int64_t pi_y = pi[y];
    int64_t sqrtz = isqrt(z);
    int64_t prime_c = nth_prime(c);

    // Find all trivial leaves: n = primes[b] * primes[l]
    // which satisfy phi(x / n), b - 1) = 1
    #pragma omp parallel for num_threads(threads) reduction(+: s2_trivial)
    for (int64_t i = 0; i < threads; i++)
    {
      int64_t start = max(prime_c, sqrtz) + 1;
      int64_t thread_interval = ceil_div(y - start, threads);
      start += thread_interval * i;
      int64_t stop = min(start + thread_interval, y);
      primesieve::iterator iter(start - 1, stop);
      T prime;

      while ((prime = iter.next_prime()) < stop)
      {
        int64_t xn = (int64_t) max(x / (prime * prime), prime);
        s2_trivial += pi_y - pi[xn];
      }
    }

    if (is_backup(get_wtime() - time))
      save_file(x, y, c, s2_trivial, time);
  }

  print("S2_trivial", s2_trivial, time);
  return s2_trivial;
}

} // namespace S2_trivial
} // namespace

namespace primecount {

int64_t S2_trivial(int64_t x,
                   int64_t y,
                   int64_t z,
                   int64_t c,
                   int threads)
{
  return S2_trivial::S2_trivial(x, y, z, c, threads);
}

#ifdef HAVE_INT128_T

int128_t S2_trivial(int128_t x,
                    int64_t y,
                    int64_t z,
                    int64_t c,
                    int threads)
{
  return S2_trivial::S2_trivial(x, y, z, c, threads);
}

#endif

} // namespace primecount
