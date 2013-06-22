#include <stdint.h>
#include <algorithm>
#include <cmath>
#include <limits>

using namespace std;

namespace {

/// This calculates the logarithmic integral using Ramanujan's fast
/// converging formula (accurate up to 10^17).
/// \mathrm{li}(x) = \gamma + \ln \ln{x} + \sum_{n=0}^{\infty} \frac{(-1)^{n-1}(\ln{x})^n}{n!2^{n-1}} \sum_{k=0}^{\lfloor{(n-1)/2}\rfloor} \frac{1}{2k+1}
// @see http://mathworld.wolfram.com/LogarithmicIntegral.html, (15)
///
long double li(long double x)
{
  long double GAMMA = 0.57721566490153286061;
  long double logx = log(x);
  long double sum = 0;
  long double inner_sum = 0;
  long double factorial = 1;
  long double p = -1;
  long double power2 = 1;

  int k = 0;
  int terms = static_cast<int>(logx * 2) + 10;

  for (int n = 1; n < terms; n++)
  {
    factorial *= n;
    p *= -logx;
    long double q = factorial * power2;
    power2 *= 2;
    for (; k <= (n - 1) / 2; k++)
      inner_sum += 1.0 / (2 * k + 1);
    sum += (p / q) * inner_sum;
  }

  return GAMMA + log(logx) + sqrt(x) * sum;
}

} // namespace

namespace primecount {

/// This calculates the offset logarithmic integral which is a very
/// accurate approximation of the number of primes below x. The
/// logarithmic integral Li(x) is larger than pi(x) below ~ 10^316.
/// \mathrm{Li}(x) = \int_2^x \frac{dt}{\log{t}}
///
int64_t Li(int64_t x)
{
  if (x < 2) return 0;
  long double n = static_cast<long double>(x);
  return static_cast<int64_t>(
      li(n) - /* li(2) = */ 1.04516);
}

/// This calculates the inverse logarithmic integral Li^{-1}(x) using
/// a binary search over the interval [n, n * log(n) * 2].
/// This function approximates the nth prime very accurately.
/// The inverse logarithmic integral Li^{-1}(x) is smaller than the
/// nth prime below ~ 10^316.
///
int64_t Li_inverse(int64_t n)
{
  double x = static_cast<double>(n);
  double logx = log(x);
  int64_t first = n;
  int64_t last = static_cast<int64_t>(x * logx * 2 + 10000);
  if (last <= first)
    last = std::numeric_limits<int64_t>::max();
  int64_t len = last - first;

  while (len != 0)
  {
    int64_t len2 = len / 2;
    int64_t guess = first + len2;
    if (n < Li(guess))
      len = len2;
    else
    {
      first = guess + 1;
      len -= len2 + 1;
    }
  }

  return first;
}

} // namespace primecount
