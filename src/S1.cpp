///
/// @file  S1.cpp
/// @brief Functions to calculate the contribution of the ordinary
///        leaves in the Lagarias-Miller-Odlyzko and Deleglise-Rivat
///        prime counting algorithms.
///
/// Copyright (C) 2015 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <S1.hpp>
#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <PhiTiny.hpp>
#include <generate.hpp>
#include <pmath.hpp>

#include <stdint.h>
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
namespace S1 {

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
               T s1,
               double time)
{
  ofstream outfile("S1.txt");

  if (outfile.is_open())
  {
    outfile << "x = " << x << endl;
    outfile << "y = " << y << endl;
    outfile << "c = " << c << endl;
    outfile << "S1 = " << s1 << endl;
    outfile << "Seconds = " << fixed << setprecision(3) << (get_wtime() - time) << endl;
    outfile << "Status = 100%" << endl;
    outfile.close();
  }
}

template <typename T>
void read_file(T x,
               int64_t y,
               int64_t c,
               T* s1,
               double* time)
{
  ifstream infile("S1.txt");

  if (infile.is_open())
  {
    try
    {
      T x2 = get_next_line(infile);
      int64_t y2 = (int64_t) get_next_line(infile);
      int64_t c2 = (int64_t) get_next_line(infile);
      int64_t s12 = (int64_t) get_next_line(infile);
      double seconds = get_next_double(infile);

      infile.close();

      if (x == x2 &&
          y == y2 &&
          c == c2)
      {
        *s1 = s12;
        *time -= seconds;

        if (print_status())
        {
          if (!print_variables())
            cout << endl;

          cout << "--- Resuming from S1.txt ---" << endl;
          cout << "x = " << x << endl;
          cout << "y = " << y << endl;
          cout << "c = " << c << endl;
          cout << "S1 = " << *s1 << endl;
          cout << "Seconds = " << seconds << endl;
          cout << "Status = 100%" << endl;
          cout << endl;
        }
      }
    }
    catch (std::exception&)
    {
      throw primecount_error("failed to read S1.txt");
    }
  }
}

/// Recursively iterate over the square free numbers coprime to the
/// first b primes and calculate the sum of the ordinary leaves.
/// This algorithm is based on section 2.2 of the paper:
/// Douglas Staple, "The Combinatorial Algorithm For Computing pi(x)",
/// arXiv:1503.01839, 6 March 2015.
///
template <int MU, typename T, typename P>
T S1(T x,
     int64_t y,
     int64_t b,
     int64_t c,
     T square_free,
     vector<P>& primes)
{
  T s1 = 0;

  for (b += 1; b < (int64_t) primes.size(); b++)
  {
    T next = square_free * primes[b];
    if (next > y) break;
    s1 += MU * phi_tiny(x / next, c);
    s1 += S1<-MU>(x, y, b, c, next, primes);
  }

  return s1;
}

/// Calculate the contribution of the ordinary leaves in parallel.
/// Run time: O(y * log(log(y))) operations.
/// Space complexity: O(y / log(y)).
///
template <typename X, typename Y>
X S1(X x,
     Y y,
     int64_t c,
     int threads)
{
  int64_t thread_threshold = ipow(10, 6);
  threads = validate_threads(threads, y, thread_threshold);
  vector<Y> primes = generate_primes<Y>(y);
  X s1 = phi_tiny(x, c);

  #pragma omp parallel for schedule(static, 1) num_threads(threads) reduction (+: s1)
  for (int64_t b = c + 1; b < (int64_t) primes.size(); b++)
  {
    s1 += -1 * phi_tiny(x / primes[b], c);
    s1 += S1<1>(x, y, b, c, (X) primes[b], primes);
  }

  return s1;
}

} // namespace S1
} // namespace

namespace primecount {

int64_t S1(int64_t x,
           int64_t y,
           int64_t c,
           int threads)
{
  print("");
  print("=== S1(x, y) ===");
  print("Computation of the ordinary leaves");
  print(x, y, c, threads);

  int64_t s1 = -1;
  double time = get_wtime();

  S1::read_file(x, y, c, &s1, &time);

  if (s1 < 0)
    s1 = S1::S1(x, y, c, threads);

  if (is_backup(get_wtime() - time))
    S1::save_file(x, y, c, s1, time);

  print("S1", s1, time);
  return s1;
}

#ifdef HAVE_INT128_T

int128_t S1(int128_t x,
            int64_t y,
            int64_t c,
            int threads)
{
  print("");
  print("=== S1(x, y) ===");
  print("Computation of the ordinary leaves");
  print(x, y, c, threads);

  double time = get_wtime();
  int128_t s1 = -1;

  S1::read_file(x, y, c, &s1, &time);

  if (s1 < 0)
  {
    // uses less memory
    if (y <= numeric_limits<uint32_t>::max())
      s1 = S1::S1(x, (uint32_t) y, c, threads);
    else
      s1 = S1::S1(x, y, c, threads);
  }

  if (is_backup((get_wtime() - time) * 1000))
    S1::save_file(x, y, c, s1, time);

  print("S1", s1, time);
  return s1;
}

#endif

} // namespace primecount
