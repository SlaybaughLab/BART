#include "utility/has_value.hpp"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

using namespace bart;

template <typename ValueType>
class UtilityHasAggregatedValueTest : public ::testing::Test {
 public:
  [[nodiscard]] auto GetValue() const -> ValueType;
};

template<>
auto UtilityHasAggregatedValueTest<double>::GetValue() const -> double {
  return test_helpers::RandomDouble(-100, 100);
}

using TestTypes = ::testing::Types<double>;
TYPED_TEST_SUITE(UtilityHasAggregatedValueTest, TestTypes);

TYPED_TEST(UtilityHasAggregatedValueTest, Set) {
  utility::HasValue<TypeParam> has_value_class;
  TypeParam value_to_set{ this->GetValue() };
  has_value_class.SetValue(value_to_set);
  EXPECT_EQ(has_value_class.value(), value_to_set);
}

TYPED_TEST(UtilityHasAggregatedValueTest, Add) {
  utility::HasValue<TypeParam> has_value_class;
  TypeParam value_to_set{ this->GetValue() }, value_to_add{ this->GetValue() };

  has_value_class.SetValue(value_to_set);
  EXPECT_EQ(has_value_class.value(), value_to_set);
  has_value_class.Add(value_to_add);
  EXPECT_EQ(has_value_class.value(), value_to_set + value_to_add);
}

} // namespace
