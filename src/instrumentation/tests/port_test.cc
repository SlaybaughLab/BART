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

TYPED_TEST(InstrumentationPortTest, Expose) {
  using InstrumentType = instrumentation::InstrumentMock<TypeParam>;
  auto instrument_ptr = std::make_shared<InstrumentType>();

  struct TestName;

  instrumentation::Port<TypeParam, TestName> port;
  port.AddInstrument(instrument_ptr);
  auto test_value = this->GetTestValue();
  EXPECT_CALL(*instrument_ptr, Read(test_value));
  port.Expose(test_value);
}

TYPED_TEST(InstrumentationPortTest, ExposeNoInstrument) {
  struct TestName;
  instrumentation::Port<TypeParam, TestName> port;
  auto test_value = this->GetTestValue();
  EXPECT_NO_THROW({
                    port.Expose(test_value);
                  });
}

TYPED_TEST(InstrumentationPortTest, GetPortHelperFunction) {
  struct PortName;
  using TestPort = instrumentation::Port<TypeParam, TestTypes>;
  using InstrumentType = instrumentation::InstrumentMock<TypeParam>;
  auto instrument_ptr = std::make_shared<InstrumentType>();

  class InstrumentedClass : public TestPort {};
  InstrumentedClass test_class;

  instrumentation::GetPort<TestPort>(test_class).AddInstrument(instrument_ptr);
  auto test_value = this->GetTestValue();
  EXPECT_CALL(*instrument_ptr, Read(test_value));
  instrumentation::GetPort<TestPort>(test_class).Expose(test_value);
}

} // namespace
