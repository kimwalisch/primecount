#include <primecount.h>
#include <primesieve/soe/PrimeSieve.h>

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <stdint.h>

namespace {

void test1()
{
  std::cout << "Testing pi(x) implementations for 0 <= x < 100000 ..." << std::endl;

  PrimeSieve ps;
  int64_t pix_primesieve = 0;
  int64_t pix_legendre = 0;
  int64_t pix_meissel = 0;
  int64_t pix_lehmer = 0;
  
  for (int64_t x = 0; x < 100000; x++)
  {
    pix_primesieve = ps.countPrimes(0, x);

    pix_legendre = primecount::pi_legendre(x);
    if (pix_legendre != pix_primesieve) {
      std::cerr << "pi_legendre(" << x << ") = "          << pix_legendre
                << " is an error, the correct result is " << pix_primesieve << std::endl;
      std::exit(1);
    }

    pix_meissel = primecount::pi_meissel(x);
    if (pix_meissel != pix_primesieve) {
      std::cerr << "pi_meissel(" << x << ") = "           << pix_meissel
                << " is an error, the correct result is " << pix_primesieve << std::endl;
      std::exit(1);
    }

    pix_lehmer = primecount::pi_lehmer(x);
    if (pix_lehmer != pix_primesieve) {
      std::cerr << "pix_lehmer(" << x << ") = "           << pix_lehmer
      << " is an error, the correct result is "           << pix_primesieve << std::endl;
      std::exit(1);
    }
  }
  std::cout << "All tests passed successfully!" << std::endl;
}

void test2()
{
  std::cout << "Randomly testing pi(x) implementations up to 2^36 ..." << std::endl;

  PrimeSieve ps;
  int64_t pix_primesieve = 0;
  int64_t pix_legendre = 0;
  int64_t pix_meissel = 0;
  int64_t pix_lehmer = 0;
  int64_t old = 0;
  srand(static_cast<unsigned int>(time(0)));

  for (int64_t x = 0; (x >> 36) == 0; old = x + 1, x += rand())
  {
    pix_primesieve += ps.countPrimes(old, x);

    pix_legendre = primecount::pi_legendre(x);
    if (pix_legendre != pix_primesieve) {
      std::cerr << "pi_legendre(" << x << ") = "          << pix_legendre
                << " is an error, the correct result is " << pix_primesieve << std::endl;
      std::exit(1);
    }

    pix_meissel = primecount::pi_meissel(x);
    if (pix_meissel != pix_primesieve) {
      std::cerr << "pi_meissel(" << x << ") = "           << pix_meissel
                << " is an error, the correct result is " << pix_primesieve << std::endl;
      std::exit(1);
    }

    pix_lehmer = primecount::pi_lehmer(x);
    if (pix_lehmer != pix_primesieve) {
      std::cerr << "pix_lehmer(" << x << ") = "           << pix_lehmer
      << " is an error, the correct result is "           << pix_primesieve << std::endl;
      std::exit(1);
    }
  }
  std::cout << "All tests passed successfully!" << std::endl;
}

} // namespace

namespace primecount {

void test()
{
  test1();
  test2();
  std::exit(0);
}

} // namespace primecount
