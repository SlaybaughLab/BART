#ifndef BART_SRC_UTILITY_HAS_VALUE_HPP_
#define BART_SRC_UTILITY_HAS_VALUE_HPP_

namespace bart::utility {

template <typename ValueType>
class HasValue {
 public:
  HasValue();
  ~HasValue() = default;
  auto Add(const ValueType& to_add) -> void;
  auto SetValue(const ValueType& to_set) -> void;
  auto value() const -> ValueType;
 private:
  ValueType value_;
};

} // namespace bart::utility

#endif //BART_SRC_UTILITY_HAS_VALUE_HPP_
