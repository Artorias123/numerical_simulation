#include "NumeralCal/Utils/Types.h"
#include <array>
#include <gtest/gtest.h>
#include <vector>

template <typename EleTy, Range T>
void checkRange1(T range, EleTy val) {
  static_assert(std::is_same_v<typename TraitRangeType<T>::type, EleTy>);
  EXPECT_EQ(range[1], val);
}

template <typename EleTy, RangeOf<EleTy> T>
void checkRange0(T range, EleTy val) {
  EXPECT_EQ(range[0], val);
}

TEST(TemplateTest, Test1) {
  std::array<float, 5> array{1, 2, 3, 4, 5};
  static_assert(std::is_same_v<decltype(*(array.begin())), float &>);
  checkRange1<float>(array, 2);

  std::vector<int> vec{2, 3, 4};
  checkRange1<int>(vec, 3);

  static_assert(RangeOf<decltype(vec), int>);
  static_assert(RangeOf<decltype(array), float>);
  checkRange0<int>(vec, 2);
  checkRange0<float>(array, 1);
}
