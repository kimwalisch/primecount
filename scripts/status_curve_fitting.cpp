// 1. Run the program: c++ status_curve_fitting.cpp -o status && ./status
// 2. Paste output into curve fitting tool: https://planetcalc.com/8735/
// 3. Select "Function must pass through particular points" with (0, 0) and (100, 100)
// 4. Select "4th order polynomial regression"
// 5. Select "Calculation precision Digits after the decimal point:" 20
// 6. Click "CALCULATE" button

#include <iostream>
#include <cmath>

void skewed_percent(double percent)
{
  double exp = 0.96;
  double base = exp + percent / (101 / (1 - exp));
  double low = std::pow(base, 100.0);
  double dividend = std::pow(base, percent) - low;
  percent = 100 - (100 * dividend / (1 - low));

  std::cout << percent << " ";
}

int main()
{
  for (double i = 0; i <= 100; i += 0.25)
    std::cout << i << " ";

  std::cout << std::endl;
  std::cout << std::endl;

  for (double i = 0; i <= 100; i += 0.25)
    skewed_percent(i);

  std::cout << std::endl;
  return 0;
}
