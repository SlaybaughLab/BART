#include "utility/has_value.hpp"

namespace bart::utility {

template <typename ValueType>
auto HasValue<ValueType>::Add(const ValueType& to_add) -> void {
  value_ += to_add;
}

template <typename ValueType>
auto HasValue<ValueType>::SetValue(const ValueType& to_set) -> void {
  value_ = to_set;
}

template <typename ValueType>
auto HasValue<ValueType>::value() const -> ValueType {
  return value_;
}
template<>
HasValue<double>::HasValue() {
  value_ = 0;
}

template class HasValue<double>;

} // namespace bart::utility