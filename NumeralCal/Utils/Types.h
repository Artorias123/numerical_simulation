#ifndef NUMERALCAL_UTILS_TYPES_H_
#define NUMERALCAL_UTILS_TYPES_H_
#include <concepts>
#include <ranges>

template <typename T>
concept Float = std::floating_point<T>;

template <typename T>
concept Range = std::ranges::range<T>;

template <typename T>
struct TraitRangeType {};

template <typename T>
concept Array =
    std::is_same_v<T, std::array<typename T::value_type, std::tuple_size_v<T>>>;

template <Range T>
struct TraitRangeType<T> {
  using type = std::remove_cvref_t<decltype(*(T{}.begin()))>;
};

template <typename T, typename U>
concept RangeOf = std::ranges::range<T> &&
    std::is_same_v<typename TraitRangeType<T>::type, U>;
#endif // NUMERALCAL_UTILS_TYPES_H_
