#include <iostream>
#include <cstdlib>

namespace primecount {

void help()
{
  std::cerr << "Usage: primecount x"                                                     << std::endl
            << "Count the primes up to x < 2^63 using Lehmer's prime counting formula." << std::endl
            << "The computational complexity is O(x/ln(x)^3) operations and"            << std::endl
            << "O(x^(1/3)/ln(x)) space."                                                << std::endl;
  std::exit(1);
}

void version()
{
  std::cout << "version 0.1" << std::endl;
  std::exit(0);
}

} // namespace primecount
