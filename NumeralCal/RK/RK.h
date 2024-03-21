#include <array>
#include <cassert>
#include <concepts>

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

template <typename T>
concept Float = std::floating_point<T>;

template <Float T, size_t N>
struct ParamA : public ParamA<T, N -1> {
  std::array<T, N> vals;
};

template <Float T>
struct ParamA<T, 1> {
  std::array<T, 1> vals;
};

namespace std {
template <size_t Idx, size_t N, Float T>
std::array<T, Idx + 1> std::get(ParamA<T, N> &v) {
  return static_cast<ParamA<T, Idx + 1>&>(v).vals;
}

template <size_t Idx, size_t N, Float T>
constexpr std::array<T, Idx + 1> std::get(const ParamA<T, N> &v) {
  return static_cast<const ParamA<T, Idx + 1>&>(v).vals;
}
} // namespace

template <Float T, size_t N>
struct RKParamUtils;

template <Float T, size_t stepN, typename... Args>
struct ButcherTableImpl;

template <Float T, size_t... Ns>
struct SequenceGenerator {
  constexpr SequenceGenerator(T, std::index_sequence<Ns...>) {}
  using tuple = ParamA<T, sizeof...(Ns)>;
  // RK 3:
  //    |
  // c2 | a21
  // c3 | a31 a32
  // ---|------------
  //    | b1  b2  b3
  // sizeof...(Ns) == 2
  // stepN = 3
  using paramType = ButcherTableImpl<T, sizeof...(Ns) + 1, std::array<T, Ns + 1>...>;
};

template <Float T, size_t N>
struct RKParamUtils {
  using tuple = typename decltype(SequenceGenerator(T{}, std::make_index_sequence<N>{}))::tuple;
  using paramType = typename decltype(SequenceGenerator(T{}, std::make_index_sequence<N - 1>{}))::paramType;
};

template <Float T, size_t stepN, typename... Args>
struct ButcherTableImpl : public RKParamUtils<T, stepN - 1>::tuple {
  using AType = RKParamUtils<T, stepN - 1>::tuple;
  using BType = std::array<T, stepN>;
  using CType = std::array<T, stepN - 1>;
  constexpr ButcherTableImpl(const BType &b,
                             const CType &c,
                             const Args&... a) : AType(a...), b(b), c(c) {}
  constexpr ButcherTableImpl(BType &&b,
                             CType &&c,
                             Args&&... a) :
    AType(std::forward<Args>(a)...),
    b(std::forward<BType>(b)),
    c(std::forward<CType>(c)) {}
  constexpr const AType &getA() const {
    return static_cast<const AType&>(*this);
  }
  constexpr const BType &getB() const {
    return b;
  }
  constexpr const CType &getC() const {
    return c;
  }

  BType b;
  CType c;
};

template <Float T, size_t stepN>
struct ButcherTable : RKParamUtils<T, stepN>::paramType {
  using Base = RKParamUtils<T, stepN>::paramType;
  using Base::Base;
  using FloatTy = T;
  inline static const size_t step = stepN;
  template <int p, int s = stepN, bool stop = (s - 2 < p)>
  constexpr bool isUseK() const {
    if constexpr (stop) {
      return Base::getB()[p] != 0;
    } else {
      // RK step is N, last A tuple index is N - 2.
      // Example:
      // RK 3 -> A: tuple<array<1>, array<2>> -> lastA = std::get<1>(A)
      return std::get<s - 2>(Base::getA())[p] != 0 || isUseK<p, s - 1>();
    }
  }
  template <int p, int s = stepN, bool stop = (s - 1 < p)>
  constexpr bool isNeedRecordK() const {
    if constexpr (stop) {
      return false;
    } else {
      // RK step is N, last A tuple index is N - 2.
      // Example:
      // RK 3 -> A: tuple<array<1>, array<2>> -> lastA = std::get<1>(A)
      return std::get<s - 2>(Base::getA())[p - 1] != 0 || isUseK<p, s - 1>();
    }
  }
};

template <typename T>
concept ButcherTableTy = std::is_same_v<T, ButcherTable<typename T::FloatTy, T::step>>;

struct ODERK {
  template <typename Fn, auto param, size_t KNum, Float T = decltype(param)::FloatTy, size_t N = decltype(param)::step>
  static T compute(T curY, T y, T x, T h, std::array<T, KNum> &k) {
    using ParamTy = decltype(param);
    static_assert(std::is_same_v<const ButcherTable<T, N>, ParamTy>);
    return run<Fn, 0, param>(curY, y, x, h, k);
  }

  template <typename Fn, size_t step, auto param, size_t skipNum = 0, Float T = decltype(param)::FloatTy, size_t N = decltype(param)::step, bool cond = (step < N)>
  constexpr static T run(T curVal, T x, T y, T h, std::array<T, N> &k) {
    using ParamTy = decltype(param);
    static_assert(std::is_same_v<const ButcherTable<T, N>, ParamTy>);
    if constexpr (cond) {
      if constexpr (param.template isUseK<step>()) {
        curVal = runStep<Fn, step, skipNum>(param.getA(), param.getB(), param.getC(), curVal, x, y, h, k);
        return run<Fn, step + 1, param, skipNum>(curVal, x, y, h, k);
      } else {
        return run<Fn, step + 1, param, skipNum + 1>(curVal, x, y, h, k);
      }
    } else {
      return curVal;
    }
  }

  template <typename Fn, ButcherTableTy ParamTy, size_t KNum, Float T = typename ParamTy::FloatTy, size_t N = ParamTy::step>
  static T compute(T curY, const ParamTy &param, T y, T x, T h, std::array<T, KNum> &k) {
    static_assert(std::is_same_v<ButcherTable<T, N>, ParamTy>);
    return run<Fn, 0>(curY, param, y, x, h, k);
  }

  template <typename Fn, size_t step, Float T, size_t N, bool cond = (step < N)>
  constexpr static T run(T curVal, const ButcherTable<T, N> &param, T x, T y, T h, std::array<T, N> &k) {
    if constexpr (cond) {
      curVal = runStep<Fn, step, 0>(param.getA(), param.getB(), param.getC(), curVal, x, y, h, k);
      return run<Fn, step + 1>(curVal, param, x, y, h, k);
    } else {
      return curVal;
    }
  }

  template <typename Fn, int step, size_t skipNum, Float T, size_t N>
  constexpr static T runStep(const ButcherTable<T, N>::AType &a,
                             const ButcherTable<T, N>::BType &b,
                             const ButcherTable<T, N>::CType &c,
                             T curVal, T x, T y, T h, std::array<T, N> &k) {
    computeK<Fn, step, skipNum>(a, c, x, y, h, k);
    return curVal + k[step - skipNum] * b[step] * h;
  }

  template <typename Fn, int step, size_t skipNum, Float T, size_t N>
  constexpr static void computeK(const ButcherTable<T, N>::AType &a,
                                 const ButcherTable<T, N>::CType &c,
                                 T x, T y, T h, std::array<T, N> &k) {
    if constexpr (step == 0) {
      k[step - skipNum] = Fn::run(x, y);
    } else {
      T curC = c[step - 1];
      k[step - skipNum] = Fn::run(x + curC * h, computeDy<0, skipNum>(y, h, std::get<step - 1>(a), k));
    }
  }

  template <int idx, size_t skipNum, Float T, size_t step, size_t N, bool cond = (idx < step)>
  constexpr static T computeDy(T y, T h, const std::array<T, step> &A, const std::array<T, N> &k) {
    if constexpr (cond) {
      return computeDy<idx + 1, skipNum>(y + A[idx] * k[idx - skipNum] * h, h, A, k);
    } else {
      return y;
    }
  }
};

struct RK : public ODERK {};
