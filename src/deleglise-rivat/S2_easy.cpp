///
/// @file  S2_easy.cpp
/// @brief Calculate the contribution of the clustered easy leaves
///        and the sparse easy leaves in parallel using OpenMP
///        (Deleglise-Rivat algorithm).
///
/// Copyright (C) 2015 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <generate.hpp>
#include <int128.hpp>
#include <min_max.hpp>
#include <pmath.hpp>
#include <S2Status.hpp>
#include <fast_div.hpp>

#include <stdint.h>
#include <algorithm>
#include <fstream>
#include <sstream>
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
namespace S2_easy {

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
               int64_t z,
               int64_t b,
               int64_t b_max,
               int64_t c,
               T s2_easy,
               double time,
               double percent)
{
  ofstream outfile("S2_easy.txt");

  if (outfile.is_open())
  {
    outfile << "x = " << x << endl;
    outfile << "y = " << y << endl;
    outfile << "z = " << z << endl;
    outfile << "b = " << b << endl;
    outfile << "b_max = " << b_max << endl;
    outfile << "c = " << c << endl;
    outfile << "S2_easy = " << s2_easy << endl;
    outfile << "Seconds = " << fixed << setprecision(3) << (get_wtime() - time) << endl;
    outfile << "Status = " << fixed << setprecision(get_status_precision(x)) << percent << '%' << endl;
    outfile.close();
  }
}

template <typename T>
void read_file(T x,
               int64_t y,
               int64_t z,
               int64_t* b,
               int64_t b_max,
               int64_t c,
               T* s2_easy,
               double* time)
{
  ifstream infile("S2_easy.txt");

  if (infile.is_open())
  {
    try
    {
      T x2 = get_next_line(infile);
      int64_t y2 = (int64_t) get_next_line(infile);
      int64_t z2 = (int64_t) get_next_line(infile);
      int64_t b2 = (int64_t) get_next_line(infile);
      int64_t b_max2 = (int64_t) get_next_line(infile);
      int64_t c2 = (int64_t) get_next_line(infile);
      T s2_easy2 = get_next_line(infile);
      double seconds = get_next_double(infile);
      double percent = get_next_double(infile);

      infile.close();

      // only resume if S2_easy.txt matches the
      // command-line values x and alpha
      if (x == x2 &&
          y == y2 &&
          z == z2 &&
          c == c2 &&
          *b < b2 &&
          b_max == b_max2)
      {
        *b = b2;
        *s2_easy = s2_easy2;
        *time -= seconds;

        if (print_status())
        {
          if (!print_variables())
            cout << endl;

          cout << "--- Resuming from S2_easy.txt ---" << endl;
          cout << "b = " << *b << endl;
          cout << "b_max = " << b_max << endl;
          cout << "S2_easy = " << *s2_easy << endl;
          cout << "Seconds = " << seconds << endl;
          cout << "Status = " << fixed << setprecision(get_status_precision(x)) << percent << '%' << endl;
          cout << endl;
        }

        // we must start at b + 1
        *b += 1;
      }
    }
    catch (std::exception&)
    {
      throw primecount_error("failed to read S2_easy.txt");
    }
  }
}

/// Calculate the contribution of the clustered easy leaves
/// and the sparse easy leaves.
/// @param T  either int64_t or uint128_t.
///
template <typename T1, typename T2>
T1 S2_easy(T1 x,
           int64_t y,
           int64_t z,
           int64_t c,
           vector<T2>& primes,
           int threads,
           double& time)
{
  T1 s2_easy = 0;
  int64_t x13 = iroot<3>(x);
  int64_t thread_threshold = 1000;
  threads = validate_threads(threads, x13, thread_threshold);

  int64_t pi_sqrty = pi_bsearch(primes, isqrt(y));
  int64_t pi_x13 = pi_bsearch(primes, x13);
  int64_t start = max(c, pi_sqrty) + 1;

  read_file(x, y, z, &start, pi_x13, c, &s2_easy, &time);

  if (start <= pi_x13)
  {
    PiTable pi(y);
    S2Status status(x);
    int64_t indexes_per_thread = 1;

    double backup_time = get_wtime();

    while (start <= pi_x13)
    {
      int64_t stop = min(start + indexes_per_thread * threads, pi_x13);

      #pragma omp parallel for schedule(dynamic, 1) num_threads(threads) reduction(+: s2_easy)
      for (int64_t b = start; b <= stop; b++)
      {
        int64_t prime = primes[b];
        T1 x2 = x / prime;
        int64_t min_trivial_leaf = min(x2 / prime, y);
        int64_t min_clustered_easy_leaf = (int64_t) isqrt(x2);
        int64_t min_sparse_easy_leaf = z / prime;
        int64_t min_hard_leaf = max(y / prime, prime);

        min_sparse_easy_leaf = max(min_sparse_easy_leaf, min_hard_leaf);
        min_clustered_easy_leaf = max(min_clustered_easy_leaf, min_hard_leaf);
        int64_t l = pi[min_trivial_leaf];
        T1 sum = 0;

        // Find all clustered easy leaves:
        // n = primes[b] * primes[l]
        // x / n <= y && phi(x / n, b - 1) == phi(x / m, b - 1)
        // where phi(x / n, b - 1) = pi(x / n) - b + 2
        while (primes[l] > min_clustered_easy_leaf)
        {
          int64_t xn = (int64_t) fast_div(x2, primes[l]);
          int64_t phi_xn = pi[xn] - b + 2;
          int64_t last_prime = primes[b + phi_xn - 1];
          int64_t xm = max((int64_t) fast_div(x2, last_prime), min_clustered_easy_leaf);
          int64_t l2 = pi[xm];
          sum += phi_xn * (l - l2);
          l = l2;
        }

        // Find all sparse easy leaves:
        // n = primes[b] * primes[l]
        // x / n <= y && phi(x / n, b - 1) = pi(x / n) - b + 2
        for (; primes[l] > min_sparse_easy_leaf; l--)
        {
          int64_t xn = (int64_t) fast_div(x2, primes[l]);
          sum += pi[xn] - b + 2;
        }

        if (print_status())
          status.print(b, pi_x13);

        s2_easy += sum;
      }

      if (is_backup(get_wtime() - backup_time))
      {
        indexes_per_thread = max(indexes_per_thread / 2, 1);
        save_file(x, y, z, stop, pi_x13, c, s2_easy, time, status.skewed_percent(stop, pi_x13));
        backup_time = get_wtime();
      }
      else
      {
        indexes_per_thread *= 2;
      }

      start = stop + 1;
    }

    if (is_backup(get_wtime() - time))
      save_file(x, y, z, pi_x13, pi_x13, c, s2_easy, time, 100);
  }

  return s2_easy;
}

} // namespace S2_easy
} // namespace

namespace primecount {

int64_t S2_easy(int64_t x,
                int64_t y,
                int64_t z,
                int64_t c,
                int threads)
{
  print("");
  print("=== S2_easy(x, y) ===");
  print("Computation of the easy special leaves");
  print(x, y, c, threads);

  double time = get_wtime();
  vector<int32_t> primes = generate_primes(y);
  int64_t s2_easy = S2_easy::S2_easy((intfast64_t) x, y, z, c, primes, threads, time);

  print("S2_easy", s2_easy, time);
  return s2_easy;
}

#ifdef HAVE_INT128_T

int128_t S2_easy(int128_t x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 int threads)
{
  print("");
  print("=== S2_easy(x, y) ===");
  print("Computation of the easy special leaves");
  print(x, y, c, threads);

  double time = get_wtime();
  int128_t s2_easy;

  // uses less memory
  if (y <= numeric_limits<uint32_t>::max())
  {
    vector<uint32_t> primes = generate_primes<uint32_t>(y);
    s2_easy = S2_easy::S2_easy((intfast128_t) x, y, z, c, primes, threads, time);
  }
  else
  {
    vector<int64_t> primes = generate_primes<int64_t>(y);
    s2_easy = S2_easy::S2_easy((intfast128_t) x, y, z, c, primes, threads, time);
  }

  print("S2_easy", s2_easy, time);
  return s2_easy;
}

#endif

} // namespace primecount
