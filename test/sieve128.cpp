#include <primecount-internal.hpp>
#include <int128_t.hpp>

#include <algorithm>
#include <iostream>

using namespace primecount;

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
    nth_prime = nth_prime_sieve_forward(n, start);
  else
    nth_prime = nth_prime_sieve_backward(-n, start);

  std::cout << "nth_prime: " << nth_prime << std::endl;

  return 0;
}
