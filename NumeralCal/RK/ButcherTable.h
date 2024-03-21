#ifndef NUMERALCAL_RK_BUTCHERTABLE_H_
#define NUMERALCAL_RK_BUTCHERTABLE_H_
#include <array>
#include <cassert>

#include "NumeralCal/Utils/Types.h"

/*
 *  R-K method:
 *  y' = f(x, y)
 *  K1 = f(x_i, y_i)
 *  K2 = f(x_i + c2 * h, y_i + a21 * K1)
 *  K3 = f(x_i + c3 * h, y_i + a31 * K1 + a32 * K2)
 *  ...
 *
 *  y_i+1 = y_i + b1 * K1 + b2 * K2 + b3 * K3 + ...
 *
 *  ButcherTable:
 *     |
 *  c2 | a21
 *  c3 | a31 a32
 *  .. | ...
 *  ---|------------
 *     | b1  b2  b3 ...
 */
namespace rk_utils {
template <Float T, size_t N, size_t Idx, bool stop>
struct ParamA;

template <Float T, size_t N, size_t Idx = 1, bool stop = (Idx == N)>
struct ParamA : public ParamA<T, N, Idx + 1> {
  using ParamTy = std::array<T, Idx>;
  template <Array... Args>
  constexpr ParamA(const ParamTy &v, const Args... args)
      : ParamA<T, N, Idx + 1>(args...), vals(v) {}
  ParamTy vals;
};

template <Float T, size_t N, size_t Idx>
struct ParamA<T, N, Idx, true> {
  using ParamTy = std::array<T, Idx>;
  constexpr ParamA(const ParamTy &v) : vals(v) {}
  ParamTy vals;
};
} // namespace rk_utils

namespace std {
template <size_t Idx, size_t N, Float T>
inline constexpr std::array<T, Idx + 1> &std::get(rk_utils::ParamA<T, N> &v) {
  return static_cast<rk_utils::ParamA<T, N, Idx + 1> &>(v).vals;
}

template <size_t Idx, size_t N, Float T>
inline constexpr const std::array<T, Idx + 1> &
std::get(const rk_utils::ParamA<T, N> &v) {
  return static_cast<const rk_utils::ParamA<T, N, Idx + 1> &>(v).vals;
}
} // namespace std

namespace rk_utils {
template <Float T, size_t N>
struct RKParamUtils;

template <Float T, size_t stepN, typename... Args>
struct ButcherTableImpl;

template <Float T, size_t... Ns>
struct RKSequenceGenerator {
  constexpr RKSequenceGenerator(T, std::index_sequence<Ns...>) {}
  using tuple = rk_utils::ParamA<T, sizeof...(Ns)>;
  // RK 3:
  //    |
  // c2 | a21
  // c3 | a31 a32
  // ---|------------
  //    | b1  b2  b3
  // sizeof...(Ns) == 2
  // stepN = 3
  using paramType =
      ButcherTableImpl<T, sizeof...(Ns) + 1, std::array<T, Ns + 1>...>;
};

template <Float T, size_t N>
struct RKParamUtils {
  using tuple = typename decltype(RKSequenceGenerator(
      T{}, std::make_index_sequence<N>{}))::tuple;
  using paramType = typename decltype(RKSequenceGenerator(
      T{}, std::make_index_sequence<N - 1>{}))::paramType;
};

template <Float T>
struct RKParamUtils<T, 0> {};

template <Float T, size_t stepN, typename... Args>
struct ButcherTableImpl : public RKParamUtils<T, stepN - 1>::tuple {
  using AType = RKParamUtils<T, stepN - 1>::tuple;
  using BType = std::array<T, stepN>;
  using CType = std::array<T, stepN - 1>;
  constexpr ButcherTableImpl(const BType &b, const CType &c, const Args &...a)
      : AType(a...), b(b), c(c) {}
  constexpr ButcherTableImpl(BType &&b, CType &&c, Args &&...a)
      : AType(std::forward<Args>(a)...), b(std::forward<BType>(b)),
        c(std::forward<CType>(c)) {}
  inline constexpr const AType &getA() const {
    return static_cast<const AType &>(*this);
  }
  inline constexpr const BType &getB() const { return b; }
  inline constexpr const CType &getC() const { return c; }

  BType b;
  CType c;
};
} // namespace rk_utils
template <Float T, size_t stepN>
struct ButcherTable;

template <typename T>
concept ButcherTableTy =
    std::is_same_v<T, ButcherTable<typename T::FloatTy, T::step>>;

template <Float T, size_t stepN>
struct ButcherTable : rk_utils::RKParamUtils<T, stepN>::paramType {
  static_assert(stepN > 1 && "ButcherTable must stepN > 1");
  using Base = rk_utils::RKParamUtils<T, stepN>::paramType;
  using Base::Base;
  using FloatTy = T;
  inline static const size_t step = stepN;
  template <int p>
  inline constexpr bool isUseK() const {
    assert(p < stepN && "p >= stepN");
    // if (a_ij == 0 with i = 0:step):
    //    no need to record Kj
    auto result = Base::getB()[p] != 0 || needRecordK<p>();
    if constexpr (p == stepN - 1) {
      return result;
    } else {
      return result || std::get<p>(Base::getA())[p] != 0;
    }
  }
  template <int p, int s = stepN>
  inline constexpr bool needRecordK() const {
    assert(p < stepN && "p >= stepN");
    if constexpr (s - 3 < p) {
      return false;
    } else {
      // RK step is N, last A tuple index is N - 2.
      // Example:
      // RK 3 -> A: tuple<array<1>, array<2>> -> lastA = std::get<1>(A)
      // if (a_ij == 0 when i - 1 != j):
      //    no need to record Kj
      return std::get<s - 2>(Base::getA())[p] != 0 || needRecordK<p, s - 1>();
    }
  }
  inline constexpr int needKNum() const { return needKNumImpl(); }

private:
  template <int p = stepN - 1>
  inline constexpr int needKNumImpl() const {
    int result{0};
    if constexpr (p > 0) {
      result = needKNumImpl<p - 1>();
    }
    if (needRecordK<p>()) {
      result += 1;
    }
    return result;
  }
};
#endif // NUMERALCAL_RK_BUTCHERTABLE_H_
