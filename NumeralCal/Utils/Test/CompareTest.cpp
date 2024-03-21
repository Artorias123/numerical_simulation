#include "NumeralCal/Utils/Compare.h"
#include <array>

#include <gtest/gtest.h>

static const std::array<float, 3> x{1, 2, 3};
static const std::array<double, 3> y{3, 5, 7};
TEST(DiffTest, Test1) {
  auto diff1 = diff<DiffType::kDiff1>(x, y);
  static_assert(std::is_same_v<decltype(diff1), double>);
  EXPECT_EQ(diff1, (2.0 + 3.0 / 2 + 4.0 / 3) / 3);
  auto diff2 = diff<DiffType::kDiff2>(x, y);
  EXPECT_EQ(diff2, std::sqrt(4.0 / 1 + 9.0 / 4 + 16.0 / 9) / 3);
  auto diff3 = diff<DiffType::kDiff3>(x, y);
  EXPECT_EQ(diff3, 2);
  auto diff4 = diff<DiffType::kDiff1, false>(x, y);
  EXPECT_EQ(diff4, 3);
  auto diff5 = diff<DiffType::kDiff2, false>(x, y);
  EXPECT_EQ(diff5, std::sqrt(29) / 3);
  auto diff6 = diff<DiffType::kDiff3, false>(x, y);
  EXPECT_EQ(diff6, 4);
}
