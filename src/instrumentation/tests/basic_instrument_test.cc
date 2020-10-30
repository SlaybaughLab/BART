#include "instrumentation/basic_instrument.h"

#include "instrumentation/outstream/tests/outstream_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

namespace instrumentation = bart::instrumentation;
using ::testing::ReturnRef;

template <typename DataType>
class InstrumentationBasicInstrumentTest : public ::testing::Test {
 public:
  using InstrumentType = instrumentation::BasicInstrument<DataType>;
  using OutstreamType = instrumentation::outstream::OutstreamMock<DataType>;
  std::unique_ptr<InstrumentType> test_instrument = nullptr;
  OutstreamType* outstream_obs_ptr_ = nullptr;
  DataType GetTestValue();
  void SetUp() override;
};

template <typename DataType>
void InstrumentationBasicInstrumentTest<DataType>::SetUp() {
  test_instrument = std::make_unique<InstrumentType>(
      std::make_unique<OutstreamType>());
  outstream_obs_ptr_ = dynamic_cast<OutstreamType*>(
      test_instrument->outstream_ptr());
}

template <>
std::string InstrumentationBasicInstrumentTest<std::string>::GetTestValue() {
  return std::string{"test string"}; }
template <>
double InstrumentationBasicInstrumentTest<double>::GetTestValue() {
  return bart::test_helpers::RandomDouble(0, 100);
}
template <>
int InstrumentationBasicInstrumentTest<int>::GetTestValue() {
  return bart::test_helpers::RandomDouble(0, 100);
}

using TestTypes = ::testing::Types<std::string, double, int>;
TYPED_TEST_SUITE(InstrumentationBasicInstrumentTest, TestTypes);

TYPED_TEST(InstrumentationBasicInstrumentTest, ConstrutorThrow) {
  EXPECT_ANY_THROW({
    instrumentation::BasicInstrument<TypeParam> test_instrument(nullptr);
  });
}

TYPED_TEST(InstrumentationBasicInstrumentTest, Getter) {
  EXPECT_NE(nullptr, this->test_instrument->outstream_ptr());
}

TYPED_TEST(InstrumentationBasicInstrumentTest, Read) {
  auto input = this->GetTestValue();
  EXPECT_CALL(*this->outstream_obs_ptr_, Output(input))
      .WillOnce(ReturnRef(*this->outstream_obs_ptr_));
  this->test_instrument->Read(input);
}


} // namespace

