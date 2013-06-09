#include "../legendre/pi.h"
#include "../meissel/pi.h"

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
  int64_t pix_legendre = 0;
  int64_t pix_meissel = 0;
  int64_t pix_primesieve = 0;

  for (int64_t x = 0; x < 100000; x++)
  {
    pix_primesieve = ps.countPrimes(0, x);

    pix_legendre = legendre::pi(x);
    if (pix_legendre != pix_primesieve) {
      std::cerr << "legendre::pi(" << x << ") = "         << pix_legendre
                << " is an error, the correct result is " << pix_primesieve << std::endl;
      std::exit(1);
    }

    pix_meissel = meissel::pi(x);
    if (pix_meissel != pix_primesieve) {
      std::cerr << "meissel::pi(" << x << ") = "          << pix_meissel
                << " is an error, the correct result is " << pix_primesieve << std::endl;
      std::exit(1);
    }
  }
  std::cout << "All tests passed successfully!" << std::endl;
}

void test2()
{
  std::cout << "Randomly testing pi(x) implementations up to 2^36 ..." << std::endl;

  PrimeSieve ps;
  int64_t pix_legendre = 0;
  int64_t pix_meissel = 0;
  int64_t pix_primesieve = 0;
  int64_t old = 0;
  srand(static_cast<unsigned int>(time(0)));

  for (int64_t x = 0; (x >> 36) == 0; old = x + 1, x += rand())
  {
    pix_primesieve += ps.countPrimes(old, x);

    pix_legendre = legendre::pi(x);
    if (pix_legendre != pix_primesieve) {
      std::cerr << "legendre::pi(" << x << ") = "         << pix_legendre
                << " is an error, the correct result is " << pix_primesieve << std::endl;
      std::exit(1);
    }

    pix_meissel = meissel::pi(x);
    if (pix_meissel != pix_primesieve) {
      std::cerr << "meissel::pi(" << x << ") = "          << pix_meissel
                << " is an error, the correct result is " << pix_primesieve << std::endl;
      std::exit(1);
    }
  }
  std::cout << "All tests passed successfully!" << std::endl;
}

} // end namespace

int main()
{
  test1();
  test2();
  return 0;
}

