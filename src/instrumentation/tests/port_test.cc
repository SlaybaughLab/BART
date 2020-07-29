#include "instrumentation/port.h"

#include "instrumentation/tests/instrument_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

namespace instrumentation = bart::instrumentation;

template <typename DataType>
class InstrumentationPortTest : public ::testing::Test {
 public:
  using InstrumentType = instrumentation::InstrumentMock<DataType>;
  DataType GetTestValue();
  void SetUp() override;
};

template <typename DataType>
void InstrumentationPortTest<DataType>::SetUp() {}


template <>
std::string InstrumentationPortTest<std::string>::GetTestValue() {
  return std::string{"test string"}; }
template <>
double InstrumentationPortTest<double>::GetTestValue() {
  return bart::test_helpers::RandomDouble(0, 100);
}
template <>
int InstrumentationPortTest<int>::GetTestValue() {
  return bart::test_helpers::RandomDouble(0, 100);
}

using TestTypes = ::testing::Types<std::string, double, int>;
TYPED_TEST_SUITE(InstrumentationPortTest, TestTypes);

TYPED_TEST(InstrumentationPortTest, AddInstrument) {
  using InstrumentType = instrumentation::InstrumentMock<TypeParam>;
  auto instrument_ptr = std::make_shared<InstrumentType>();

  struct TestName;

  instrumentation::Port<TypeParam, TestName> port;
  port.AddInstrument(instrument_ptr);
  EXPECT_NE(port.instrument_ptr(), nullptr);
}

} // namespace
