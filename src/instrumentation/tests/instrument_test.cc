#include "instrumentation/instrument.h"

#include <string>
#include <utility>

#include "instrumentation/converter/tests/converter_mock.h"
#include "instrumentation/outstream/tests/outstream_mock.h"
#include "convergence/status.hpp"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

using namespace bart;

using ::testing::Return, ::testing::ReturnRef;

template <typename TypePair>
class InstrumentationInstrumentTest : public ::testing::Test {
 public:
  using InputType = typename TypePair::first_type;
  using OutputType = typename TypePair::second_type;
  using ConverterType = typename instrumentation::converter::ConverterMock<InputType, OutputType>;
  using OutstreamType = typename instrumentation::outstream::OutstreamMock<OutputType>;
  using InstrumentType = typename instrumentation::Instrument<InputType, OutputType>;

  std::unique_ptr<InstrumentType> test_instrument = nullptr;

  ConverterType* converter_obs_ptr_= nullptr;
  OutstreamType* outstream_obs_ptr_ = nullptr;

  InputType GetInput() const;
  OutputType GetOutput() const;

  void SetUp() override;
};

template <typename TypePair>
void InstrumentationInstrumentTest<TypePair>::SetUp() {
  auto converter_ptr = std::make_unique<ConverterType>();
  converter_obs_ptr_ = converter_ptr.get();
  auto outstream_ptr = std::make_unique<OutstreamType>();
  outstream_obs_ptr_ = outstream_ptr.get();

  test_instrument = std::make_unique<InstrumentType>(std::move(converter_ptr),
                                                     std::move(outstream_ptr));
}

template <>
convergence::Status InstrumentationInstrumentTest<std::pair<convergence::Status, std::string>>::GetInput() const {
  convergence::Status test_status;
  test_status.is_complete = false;
  test_status.delta = test_helpers::RandomDouble(1e-10, 1e-6);
  test_status.failed_index = test_helpers::RandomDouble(0, 10);
  test_status.iteration_number = test_helpers::RandomDouble(0, 1000);
  test_status.max_iterations = test_helpers::RandomDouble(1000, 10000);
  return test_status;
}

template <>
std::string InstrumentationInstrumentTest<std::pair<convergence::Status, std::string>>::GetOutput() const {
  return std::string{"test_string"};
}

using PairTypes = ::testing::Types<
    std::pair<convergence::Status, std::string>
>;

TYPED_TEST_SUITE(InstrumentationInstrumentTest, PairTypes);

TYPED_TEST(InstrumentationInstrumentTest, DepdendencyGetters) {
  EXPECT_NE(nullptr, this->test_instrument->converter_ptr());
  EXPECT_NE(nullptr, this->test_instrument->outstream_ptr());
}

TYPED_TEST(InstrumentationInstrumentTest, NullDepdendencies) {
  using InputType = typename TypeParam::first_type;
  using OutputType = typename TypeParam::second_type;
  using ConverterType = typename instrumentation::converter::ConverterMock<InputType, OutputType>;
  using OutstreamType = typename instrumentation::outstream::OutstreamMock<OutputType>;
  using InstrumentType = typename instrumentation::Instrument<InputType, OutputType>;

  EXPECT_ANY_THROW({
    InstrumentType test_instrument(std::make_unique<ConverterType>(), nullptr);
  });
  EXPECT_ANY_THROW({
    InstrumentType test_instrument(nullptr, std::make_unique<OutstreamType>());
  });
}

TYPED_TEST(InstrumentationInstrumentTest, ReadTest) {
  const auto input = this->GetInput();
  const auto output = this->GetOutput();

  EXPECT_CALL(*this->converter_obs_ptr_, Convert(input))
      .WillOnce(Return(output));
  EXPECT_CALL(*this->outstream_obs_ptr_, Output(output))
      .WillOnce(ReturnRef(*this->outstream_obs_ptr_));

  this->test_instrument->Read(input);
}

} // namespace
