#include "utility/to_string.hpp"

#include "domain/mesh/factory.hpp"
#include "instrumentation/converter/factory.hpp"
#include "instrumentation/outstream/factory.hpp"
#include "solver/group/factory.hpp"
#include "solver/linear/factory.hpp"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

namespace utility = bart::utility;
namespace converter = bart::instrumentation::converter;
namespace test_helpers = bart::test_helpers;

using MeshName = bart::domain::mesh::MeshName;
using ConverterName = converter::ConverterName;
using OutstreamName = bart::instrumentation::outstream::OutstreamName;
using GroupSolverName = bart::solver::group::GroupSolverName;
using LinearSolverName = bart::solver::linear::LinearSolverName;

struct NoStringConversionType{};

template <typename T>
class UtilityToStringTest : public ::testing::Test {
 public:
  [[nodiscard]] auto GetValue() const -> T;
};

using TestTypes = ::testing::Types<int, double, ConverterName, MeshName, OutstreamName, GroupSolverName, LinearSolverName, NoStringConversionType>;
TYPED_TEST_SUITE(UtilityToStringTest, TestTypes);

template <>
auto UtilityToStringTest<int>::GetValue() const -> int { return test_helpers::RandomInt(-100, 100); }

template <>
auto UtilityToStringTest<NoStringConversionType>::GetValue() const -> NoStringConversionType {
  return NoStringConversionType(); }

template <>
auto UtilityToStringTest<double>::GetValue() const -> double { return test_helpers::RandomDouble(-100, 100); }

template <>
auto UtilityToStringTest<ConverterName>::GetValue() const -> ConverterName {
  return static_cast<ConverterName>(test_helpers::RandomInt(0, 10)); }

template <> auto UtilityToStringTest<MeshName>::GetValue() const -> MeshName { return MeshName::kCartesian; }

template <>
auto UtilityToStringTest<OutstreamName>::GetValue() const -> OutstreamName {
  return OutstreamName::kToConditionalOstream; }

template <>
auto UtilityToStringTest<GroupSolverName>::GetValue() const -> GroupSolverName {
  return GroupSolverName::kDefaultImplementation; }

template <>
auto UtilityToStringTest<LinearSolverName>::GetValue() const -> LinearSolverName {
  return LinearSolverName::kGMRES; }

TYPED_TEST(UtilityToStringTest, ToStringReturnsSomething) {
  EXPECT_NE(utility::to_string(this->GetValue()).size(), 0);
}

} // namespace
