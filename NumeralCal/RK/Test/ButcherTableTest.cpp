#include "NumeralCal/RK/ButcherTable.h"

#include "gtest/gtest.h"

TEST(TestButcherTable, StaticTest) {
  //  constexpr ButcherTable<float, 1> rk1 = {{}, {2}, {}}
  constexpr ButcherTable<double, 2> rk2 = {{1, 1}, {1}, {1}};
  static_assert(rk2.isUseK<0>());
  static_assert(rk2.isUseK<1>());
  //  static_assert(rk2.isUseK<2>());
  static_assert(!rk2.needRecordK<0>());
  static_assert(!rk2.needRecordK<1>());
  //  static_assert(rk2.needRecordK<1>());
  constexpr ButcherTable<double, 3> rk3 = {{1, 1, 1}, {1, 1}, {1}, {1, 1}};
  static_assert(rk3.isUseK<0>());
  static_assert(rk3.isUseK<1>());
  static_assert(rk3.isUseK<2>());
  static_assert(rk3.needRecordK<0>());
  static_assert(!rk3.needRecordK<1>());
  static_assert(!rk3.needRecordK<2>());

  constexpr ButcherTable<double, 2> rk2_1 = {{0, 0}, {1}, {0}};
  static_assert(!rk2_1.isUseK<0>());
  static_assert(!rk2_1.isUseK<1>());
  static_assert(!rk2_1.needRecordK<0>());
  static_assert(!rk2_1.needRecordK<1>());

  constexpr ButcherTable<double, 3> rk3_1 = {{0, 0, 0}, {0, 0}, {0}, {0, 0}};
  static_assert(!rk3_1.isUseK<0>());
  static_assert(!rk3_1.isUseK<1>());
  static_assert(!rk3_1.isUseK<2>());
  static_assert(!rk3_1.needRecordK<0>());
  static_assert(!rk3_1.needRecordK<1>());
  static_assert(!rk3_1.needRecordK<2>());

  constexpr ButcherTable<float, 4> rk4 = {
      {1, 1, 1, 1}, {1, 1, 1}, {1}, {0, 1}, {0, 0, 1}};
  static_assert(rk4.isUseK<0>());
  static_assert(rk4.isUseK<1>());
  static_assert(rk4.isUseK<2>());
  static_assert(rk4.isUseK<3>());
  static_assert(!rk4.needRecordK<0>());
  static_assert(!rk4.needRecordK<1>());
  static_assert(!rk4.needRecordK<2>());
  static_assert(!rk4.needRecordK<3>());
  static_assert(rk4.needKNum() == 0);

  static_assert(sizeof(ButcherTable<float, 2>) == 16);
  static_assert(sizeof(ButcherTable<float, 4>) == 52);
  static_assert(sizeof(ButcherTable<double, 4>) == 104);
}

TEST(TestButcherTable, DynamicTest) {
  constexpr ButcherTable<double, 2> rk2 = {{1, 1}, {1}, {1}};
  EXPECT_TRUE(rk2.isUseK<0>());
  EXPECT_TRUE(rk2.isUseK<1>());
  EXPECT_TRUE(!rk2.needRecordK<0>());
  EXPECT_TRUE(!rk2.needRecordK<1>());

  constexpr ButcherTable<double, 3> rk3 = {{1, 1, 1}, {1, 1}, {1}, {1, 1}};
  EXPECT_TRUE(rk3.isUseK<0>());
  EXPECT_TRUE(rk3.isUseK<1>());
  EXPECT_TRUE(rk3.isUseK<2>());
  EXPECT_TRUE(rk3.needRecordK<0>());
  EXPECT_TRUE(!rk3.needRecordK<1>());
  EXPECT_TRUE(!rk3.needRecordK<2>());

  constexpr ButcherTable<double, 2> rk2_1 = {{0, 0}, {1}, {0}};
  EXPECT_TRUE(!rk2_1.isUseK<0>());
  EXPECT_TRUE(!rk2_1.isUseK<1>());
  EXPECT_TRUE(!rk2_1.needRecordK<0>());
  EXPECT_TRUE(!rk2_1.needRecordK<1>());

  ButcherTable<double, 3> rk3_1 = {{0, 0, 0}, {0, 0}, {0}, {0, 0}};
  EXPECT_TRUE(!rk3_1.isUseK<0>());
  EXPECT_TRUE(!rk3_1.isUseK<1>());
  EXPECT_TRUE(!rk3_1.isUseK<2>());
  EXPECT_TRUE(!rk3_1.needRecordK<0>());
  EXPECT_TRUE(!rk3_1.needRecordK<1>());
  EXPECT_TRUE(!rk3_1.needRecordK<2>());

  ButcherTable<float, 4> rk4 = {
      {1, 1, 1, 1}, {1, 1, 1}, {1}, {0, 1}, {0, 0, 1}};
  EXPECT_TRUE(rk4.isUseK<0>());
  EXPECT_TRUE(rk4.isUseK<1>());
  EXPECT_TRUE(rk4.isUseK<2>());
  EXPECT_TRUE(rk4.isUseK<3>());
  EXPECT_TRUE(!rk4.needRecordK<0>());
  EXPECT_TRUE(!rk4.needRecordK<1>());
  EXPECT_TRUE(!rk4.needRecordK<2>());
  EXPECT_TRUE(!rk4.needRecordK<3>());
}
