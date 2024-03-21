#ifndef NUMERALCAL_UTILS_COMPARE_H_
#define NUMERALCAL_UTILS_COMPARE_H_
#include <cassert>
#include <cmath>
#include <concepts>
#include <ranges>

#include "NumeralCal/Utils/Types.h"

enum class DiffType {
  kDiff1 = 1,
  kDiff2 = 2,
  kDiff3 = 3,
};

struct DiffUtils {
  template <DiffType type, Float DT, Float T>
  static inline DT accumulate(DT diff, T value) {
    if constexpr (type == DiffType::kDiff1 || type == DiffType::kDiff2) {
      return diff + value;
    } else if constexpr (type == DiffType::kDiff3) {
      return std::max(diff, value);
    }
  }

  template <DiffType type, Float RT>
  static inline RT postProcess(RT d, size_t size) {
    if constexpr (type == DiffType::kDiff1) {
      return d / size;
    } else if constexpr (type == DiffType::kDiff2) {
      return std::sqrt(d) / size;
    } else {
      return d;
    }
  }
};

template <DiffType type, bool relative, Float T1, Float T2,
          Float RT = decltype(T1{} + T2{})>
inline RT diff(T1 expected, T2 result) {
  RT d;
  if constexpr (type == DiffType::kDiff1 || type == DiffType::kDiff3) {
    d = std::abs(expected - result);
    if constexpr (relative) {
      d = d / std::abs(expected);
    }
  } else if constexpr (type == DiffType::kDiff2) {
    d = expected - result;
    d = d * d;
    if constexpr (relative) {
      d = d / (expected * expected);
    }
  }
  return d;
}

template <DiffType type, bool relative = true, Range T1, Range T2,
          Float RT = decltype(typename TraitRangeType<T1>::type{} +
                              typename TraitRangeType<T2>::type{})>
inline RT diff(T1 expected, T2 result) {
  using EleTy1 = typename TraitRangeType<T1>::type;
  using EleTy2 = typename TraitRangeType<T2>::type;
  static_assert(Float<EleTy1> && Float<EleTy2> &&
                std::is_same_v<RT, decltype(EleTy1{} + EleTy2{})>);
  assert(std::ranges::size(expected) == std::ranges::size(result));
  RT resultDiff{0};
  int64_t size = std::ranges::ssize(expected);
  for (auto idx : std::views::iota(0, size)) {
    auto curDiff = diff<type, relative>(expected[idx], result[idx]);
    resultDiff = DiffUtils::accumulate<type>(resultDiff, curDiff);
  }
  return DiffUtils::postProcess<type>(resultDiff, size);
}
#endif // NUMERALCAL_UTILS_COMPARE_H_
