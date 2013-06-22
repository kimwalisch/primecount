#include <stdint.h>
#include <algorithm>
#include <cmath>

using namespace std;

namespace {

/// This calculates the logarithmic integral using Ramanujan's fast
/// converging formula (accurate up to 10^17).
/// \mathrm{li}(x) = \gamma + \ln \ln{x} + \sum_{n=0}^{\infty} \frac{(-1)^{n-1}(\ln{x})^n}{n!2^{n-1}} \sum_{k=0}^{\lfloor{(n-1)/2}\rfloor} \frac{1}{2k+1}
// @see http://mathworld.wolfram.com/LogarithmicIntegral.html, (15)
///
double li(double x)
{
  long double GAMMA = 0.57721566490153286061;
  long double z = x;
  long double logz = log(z);
  long double sum = 0;
  long double inner_sum = 0;
  long double factorial = 1;
  long double p = -1;
  long double power2 = 1;

  int k = 0;
  int terms = static_cast<int>(logz * 2) + 10;

  for (int n = 1; n < terms; n++)
  {
    factorial *= n;
    p *= -logz;
    long double q = factorial * power2;
    power2 *= 2;
    for (; k <= (n - 1) / 2; k++)
      inner_sum += 1.0 / (2 * k + 1);
    sum += (p / q) * inner_sum;
  }

  return static_cast<double>(GAMMA + log(logz) + sqrt(z) * sum);
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
  return static_cast<int64_t>(li(x) - /* li(2) = */ 1.04516);
}

} // namespace primecount
