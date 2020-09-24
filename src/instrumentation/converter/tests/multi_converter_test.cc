#include "instrumentation/converter/multi_converter.h"

#include <memory>

#include "instrumentation/converter/tests/converter_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

namespace instrumentation = bart::instrumentation;
namespace test_helpers = bart::test_helpers;
using ::testing::Return;

class InstrumentationConverterMultiConverterTest : public ::testing::Test {
 public:
  using InputType = int;
  using IntermediateType = std::string;
  using OutputType = std::vector<double>;
  using MultiConverter = instrumentation::converter::MultiConverter<InputType, IntermediateType, OutputType>;
  using FirstStageConverter = instrumentation::converter::ConverterMock<InputType, IntermediateType>;
  using SecondStageConverter = instrumentation::converter::ConverterMock<IntermediateType, OutputType>;

  std::unique_ptr<MultiConverter> test_converter;

  FirstStageConverter* first_stage_converter_obs_ptr_;
  SecondStageConverter* second_stage_converter_obs_ptr_;
  void SetUp() override;
};

void InstrumentationConverterMultiConverterTest::SetUp() {
  auto first_stage_converter_ptr = std::make_unique<FirstStageConverter>();
  first_stage_converter_obs_ptr_ = first_stage_converter_ptr.get();
  auto second_stage_converter_ptr = std::make_unique<SecondStageConverter>();
  second_stage_converter_obs_ptr_ = second_stage_converter_ptr.get();
  test_converter = std::make_unique<MultiConverter>(
      std::move(first_stage_converter_ptr),
      std::move(second_stage_converter_ptr));
}

TEST_F(InstrumentationConverterMultiConverterTest, DependencyNullptrs) {
  const int n_dependencies{2};
  for (int i = 0; i < n_dependencies; ++i) {
      auto first_stage_converter_ptr = std::make_unique<FirstStageConverter>();
      auto second_stage_converter_ptr = std::make_unique<SecondStageConverter>();
      switch (i) {
        case 1: {
          second_stage_converter_ptr = nullptr;
          break;
        }
        default:
        case 0: {
          first_stage_converter_ptr = nullptr;
          break;
        }
      }
      EXPECT_ANY_THROW({
        MultiConverter test_converter(std::move(first_stage_converter_ptr),
                                      std::move(second_stage_converter_ptr));
      });
  }
}

TEST_F(InstrumentationConverterMultiConverterTest, DependencyGetters) {
  EXPECT_NE(nullptr, test_converter->second_stage_converter_ptr());
  EXPECT_NE(nullptr, test_converter->first_stage_converter_ptr());
}

TEST_F(InstrumentationConverterMultiConverterTest, Convert) {
  const int input{test_helpers::RandomInt(-100, 100)};
  const std::string intermediate{"intermediate_string"};
  const std::vector<double> output{
    test_helpers::RandomVector(test_helpers::RandomInt(5, 10), -100, 100)};
  EXPECT_CALL(*first_stage_converter_obs_ptr_, Convert(input))
      .WillOnce(Return(intermediate));
  EXPECT_CALL(*second_stage_converter_obs_ptr_, Convert(intermediate))
      .WillOnce(Return(output));
  auto test_output = test_converter->Convert(input);
  EXPECT_EQ(test_output, output);
}





} // namespace