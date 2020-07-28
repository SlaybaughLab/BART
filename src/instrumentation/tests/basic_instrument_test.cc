#include "instrumentation/basic_instrument.h"

#include "instrumentation/output/tests/output_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

namespace instrumentation = bart::instrumentation;

template <typename DataType>
class InstrumentationBasicInstrumentTest : public ::testing::Test {
 public:
  using InstrumentType = instrumentation::BasicInstrument<DataType>;
  using OutputterType = instrumentation::output::OutputMock<DataType>;
  std::unique_ptr<InstrumentType> test_instrument = nullptr;
  void SetUp() override;
};

template <typename DataType>
void InstrumentationBasicInstrumentTest<DataType>::SetUp() {
  test_instrument = std::make_unique<InstrumentType>(
      std::make_unique<OutputterType>());
}

using TestTypes = ::testing::Types<std::string, double, int>;
TYPED_TEST_SUITE(InstrumentationBasicInstrumentTest, TestTypes);

TYPED_TEST(InstrumentationBasicInstrumentTest, ConstrutorThrow) {
  EXPECT_ANY_THROW({
    instrumentation::BasicInstrument<TypeParam> test_instrument(nullptr);
  });
}

TYPED_TEST(InstrumentationBasicInstrumentTest, Getter) {
  EXPECT_NE(nullptr, this->test_instrument->outputter_ptr());
}


} // namespace

