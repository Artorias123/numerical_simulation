#include <iostream>
#include "NumeralCal/RK/RK.h"

#include "gtest/gtest.h"
template <Float T>
struct F {
static T run(T x, T y) {
  return y;
}
};

TEST(TestRK, Dy_y) {
  float x = 0, y = 1, h = 0.01;
  std::array<float, 2> k;
  constexpr ButcherTable<float, 2> rk2 = {{0.5, 0.5}, {1}, {1}};
  std::cout << x << ", " << y << std::endl;

  for (int i = 0; i < 100; i++) {
    y = RK::compute<F<float>, rk2>(y, x, y, h, k);
    x += h;
    std::cout << x << ", " << y << std::endl;
  }

  for (int i = 0; i < 100; i++) {
    y = RK::compute<F<float>>(y, rk2, x, y, h, k);
    x += h;
    std::cout << x << ", " << y << std::endl;
  }
}
