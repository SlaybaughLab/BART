#include "instrumentation/converter/convergence_to_string.h"

#include <iomanip>
#include <sstream>

#include "convergence/status.h"
#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class InstrumentationConverterConvergenceToStringTest : public ::testing::Test {
 public:
  using ConverterType = instrumentation::converter::ConvergenceToString;
  convergence::Status test_status;
  std::string GetExpectedOutput(convergence::Status&, ConverterType&) const;
  void SetUp() override;
};

void InstrumentationConverterConvergenceToStringTest::SetUp() {
  test_status.is_complete = false;
  test_status.iteration_number = test_helpers::RandomDouble(0, 1000);
  test_status.max_iterations = test_helpers::RandomDouble(1000, 10000);
  test_status.delta = test_helpers::RandomDouble(1e-10, 1e-6);
  test_status.failed_index = test_helpers::RandomDouble(0, 10);
}

std::string InstrumentationConverterConvergenceToStringTest::GetExpectedOutput(
    convergence::Status& status, ConverterType& converter) const {
  using OutputTerm = ConverterType::OutputTerm;
  const auto output_term_to_string_map = converter.output_term_to_string_map();
  auto output = converter.output_format();
  // Iteration number
  std::string iteration_num_string = output_term_to_string_map.at(OutputTerm::kIterationNum);
  output.replace(output.find(iteration_num_string),
                 iteration_num_string.size(),
                 std::to_string(test_status.iteration_number));
  // Max iterations
  std::string iteration_max_string = output_term_to_string_map.at(OutputTerm::kIterationMax);
  output.replace(output.find(iteration_max_string),
                 iteration_max_string.size(),
                 std::to_string(test_status.max_iterations));
  // Delta
  std::ostringstream delta_stream;
  delta_stream << test_status.delta.value();
  std::string delta_string = output_term_to_string_map.at(OutputTerm::kDelta);
  output.replace(output.find(delta_string),
                 delta_string.size(),
                 delta_stream.str());
  // Index
  std::string index_string = output_term_to_string_map.at(OutputTerm::kIndex);
  output.replace(output.find(index_string),
                 index_string.size(),
                 std::to_string(test_status.failed_index.value()));

  return output;
}

TEST_F(InstrumentationConverterConvergenceToStringTest, Constructor) {
  using ConverterType = instrumentation::converter::ConvergenceToString;
  EXPECT_NO_THROW({
    ConverterType test_converter;
  });
}

TEST_F(InstrumentationConverterConvergenceToStringTest, DefaultString) {
  ConverterType test_converter;
  EXPECT_EQ(test_converter.Convert(test_status),
            GetExpectedOutput(test_status, test_converter));
}

} // namespace
