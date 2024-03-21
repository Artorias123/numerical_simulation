#ifndef NUMERALCAL_RK_RK_H_
#define NUMERALCAL_RK_RK_H_
#include <array>
#include <cassert>

#include "NumeralCal/RK/ButcherTable.h"
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
template <Float T, size_t stepN>
struct ButcherTable;
struct DynamicRK;
struct StaticRK;

struct RKCore {
  template <ButcherTableTy T>
  inline static constexpr int needKArraySize(const T &table) {
    return std::max(1, table.needKNum());
  }
  template <typename Fn, int step, Float T, size_t N>
  inline static constexpr void runStep(const ButcherTable<T, N>::AType &a,
                                       const ButcherTable<T, N>::BType &b,
                                       const ButcherTable<T, N>::CType &c, T x,
                                       T y, T h, T &curVal,
                                       std::array<T, N> &k) {
    if constexpr (step == 0) {
      k[step] = Fn::compute(x, y);
    } else {
      T curC = c[step - 1];
      k[step] = Fn::compute(x + curC * h,
                            computeDy<0>(y, h, std::get<step - 1>(a), k));
    }
    curVal += k[step] * b[step] * h;
  }

  template <int idx, Float T, size_t step, size_t N>
  inline static constexpr T computeDy(T y, T h, const std::array<T, step> &A,
                                      const std::array<T, N> &k) {
    if constexpr (idx < step) {
      return computeDy<idx + 1>(y + A[idx] * k[idx] * h, h, A, k);
    } else {
      return y;
    }
  }
};

struct DynamicRK : private RKCore {
  template <typename Fn, ButcherTableTy ParamTy, size_t KNum,
            Float T = typename ParamTy::FloatTy, size_t N = ParamTy::step>
  inline static constexpr T compute(T curY, const ParamTy &param, T y, T x, T h,
                                    std::array<T, KNum> &k) {
    static_assert(std::is_same_v<ButcherTable<T, N>, ParamTy>);
    return run<Fn, 0>(curY, param, y, x, h, k);
  }

  template <typename Fn, size_t step, Float T, size_t N>
  inline static constexpr T run(T curVal, const ButcherTable<T, N> &param, T x,
                                T y, T h, std::array<T, N> &k) {
    if constexpr (step < N) {
      runStep<Fn, step>(param.getA(), param.getB(), param.getC(), x, y, h,
                        curVal, k);
      return run<Fn, step + 1>(curVal, param, x, y, h, k);
    } else {
      return curVal;
    }
  }
};

struct StaticRK : private RKCore {
  template <typename Fn, auto param, size_t KNum,
            Float T = decltype(param)::FloatTy,
            size_t N = decltype(param)::step>
  inline static constexpr T compute(T curY, T y, T x, T h,
                                    std::array<T, KNum> &k) {
    using ParamTy = decltype(param);
    static_assert(std::is_same_v<const ButcherTable<T, N>, ParamTy>);
    return run<Fn, 0, param>(curY, y, x, h, k);
  }

  template <typename Fn, size_t step, auto param,
            Float T = decltype(param)::FloatTy,
            size_t N = decltype(param)::step>
  inline static constexpr T run(T curVal, T x, T y, T h, std::array<T, N> &k) {
    using ParamTy = decltype(param);
    static_assert(std::is_same_v<const ButcherTable<T, N>, ParamTy>);
    if constexpr (step < N) {
      if constexpr (param.template isUseK<step>()) {
        runStep<Fn, step>(param.getA(), param.getB(), param.getC(), x, y, h,
                          curVal, k);
        return run<Fn, step + 1, param>(curVal, x, y, h, k);
      } else {
        return run<Fn, step + 1, param>(curVal, x, y, h, k);
      }
    } else {
      return curVal;
    }
  }
};
#endif // NUMERALCAL_RK_RK_H_
