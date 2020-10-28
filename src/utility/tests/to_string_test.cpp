#include "utility/to_string.hpp"

#include "domain/mesh/factory.hpp"
#include "instrumentation/converter/factory.hpp"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

namespace utility = bart::utility;
namespace converter = bart::instrumentation::converter;
namespace test_helpers = bart::test_helpers;

using MeshName = bart::domain::mesh::MeshName;
using ConverterName = converter::ConverterName;

template <typename T>
class UtilityToStringTest : public ::testing::Test {
 public:
  [[nodiscard]] auto GetValue() const -> T;
};

template <>
auto UtilityToStringTest<int>::GetValue() const -> int {
  return test_helpers::RandomInt(-100, 100);
}

template <>
auto UtilityToStringTest<double>::GetValue() const -> double {
  return test_helpers::RandomDouble(-100, 100);
}

template <>
auto UtilityToStringTest<ConverterName>::GetValue() const -> ConverterName {
  return static_cast<ConverterName>(test_helpers::RandomInt(0, 10));
}

template <>
auto UtilityToStringTest<MeshName>::GetValue() const -> MeshName {
  return MeshName::kCartesian;
}

using TestTypes = ::testing::Types<int, double, ConverterName, MeshName>;
TYPED_TEST_SUITE(UtilityToStringTest, TestTypes);

TYPED_TEST(UtilityToStringTest, ToStringReturnsSomething) {
  EXPECT_NE(utility::to_string(this->GetValue()).size(), 0);
}

} // namespace
