#include "instrumentation/instrument.h"

#include <string>
#include <utility>

#include "instrumentation/converter/tests/converter_mock.h"
#include "instrumentation/output/tests/output_mock.h"
#include "convergence/status.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename TypePair>
class InstrumentationInstrumentTest : public ::testing::Test {
 public:
  using InputType = typename TypePair::first_type;
  using OutputType = typename TypePair::second_type;
  using ConverterType = typename instrumentation::converter::ConverterMock<InputType, OutputType>;
  using OutputterType = typename instrumentation::output::OutputMock<OutputType>;
  using InstrumentType = typename instrumentation::Instrument<InputType, OutputType>;

  std::unique_ptr<InstrumentType> test_instrument = nullptr;

  ConverterType* converter_obs_ptr_= nullptr;
  OutputterType* outputter_obs_ptr_ = nullptr;
  void SetUp() override;
};

template <typename TypePair>
void InstrumentationInstrumentTest<TypePair>::SetUp() {
  auto converter_ptr = std::make_unique<ConverterType>();
  converter_obs_ptr_ = converter_ptr.get();
  auto outputter_ptr = std::make_unique<OutputterType>();
  outputter_obs_ptr_ = outputter_ptr.get();

  test_instrument = std::make_unique<InstrumentType>(std::move(converter_ptr),
                                                     std::move(outputter_ptr));
}

using PairTypes = ::testing::Types<
    std::pair<convergence::Status, std::string>
>;

TYPED_TEST_SUITE(InstrumentationInstrumentTest, PairTypes);

TYPED_TEST(InstrumentationInstrumentTest, DepdendencyGetters) {
  EXPECT_NE(nullptr, this->test_instrument->converter_ptr());
  EXPECT_NE(nullptr, this->test_instrument->outputter_ptr());
}

} // namespace
