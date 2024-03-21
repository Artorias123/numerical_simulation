#include "NumeralCal/RK/RK.h"
#include "NumeralCal/Utils/Compare.h"
#include <cmath>

#include "gtest/gtest.h"
template <Float T>
struct F {
  static T compute(T x, T y) { return y; }
};

template <int n>
constexpr std::array<float, n> expVec(float x0, float dx) {
  std::array<float, n> result;
  for (auto i = 0; i < n; i++) {
    result[i] = std::exp(x0);
    x0 += dx;
  }
  return result;
}

TEST(TestRK, CoreTest) {
  constexpr ButcherTable<float, 4> rk4 = {
      {1, 1, 1, 1}, {1, 1, 1}, {1}, {0, 1}, {0, 0, 1}};
  static_assert(RKCore::needKArraySize(rk4) == 1);
}

TEST(TestRK, Dy_y) {
  float x = 0, y = 1;
  constexpr float h = 0.01;
  std::array<float, 2> k;
  constexpr ButcherTable<float, 2> rk2 = {{0.5, 0.5}, {1}, {1}};
  constexpr auto expected0 = expVec<100>(0, h);

  std::array<float, 100> result;
  for (int i = 0; i < 100; i++) {
    result[i] = y;
    y = StaticRK::compute<F<float>, rk2>(y, x, y, h, k);
    x += h;
  }
  EXPECT_LE((diff<DiffType::kDiff3, false>(expected0, result)), 5e-5);

  constexpr auto expected1 = expVec<100>(1, h);
  for (int i = 0; i < 100; i++) {
    result[i] = y;
    y = DynamicRK::compute<F<float>>(y, rk2, x, y, h, k);
    x += h;
  }
  EXPECT_LE((diff<DiffType::kDiff3, false>(expected1, result)), 5e-4);
}
