#include "instrumentation/output/to_conditional_ostream.h"

#include <deal.II/base/conditional_ostream.h>
#include <memory>
#include <sstream>

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class InstrumentationOutputToConditionalOstreamTest : public ::testing::Test {
 public:
  using OutputType = instrumentation::output::ToConditionalOstream;
  using ConditionalOstreamType = dealii::ConditionalOStream;
  std::unique_ptr<OutputType> test_outputter;
  std::ostringstream output_stream;
  void SetUp() override;
};

void InstrumentationOutputToConditionalOstreamTest::SetUp() {
  test_outputter = std::make_unique<OutputType>(
      std::move(std::make_unique<ConditionalOstreamType>(output_stream, true)));
}

TEST_F(InstrumentationOutputToConditionalOstreamTest, Constructor) {
  std::ostringstream string_stream;
  EXPECT_NO_THROW({
    instrumentation::output::ToConditionalOstream test(std::move(
        std::make_unique<dealii::ConditionalOStream>(string_stream, true)));
    });
}

TEST_F(InstrumentationOutputToConditionalOstreamTest, StreamGetter) {
  EXPECT_NE(nullptr, test_outputter->conditional_ostream_ptr());
}

TEST_F(InstrumentationOutputToConditionalOstreamTest, OutputString) {
  std::string test_string{"this is a test string\n"},
      second_string{"so is this"};
  test_outputter->Output(test_string).Output(second_string);
  EXPECT_EQ(test_string + second_string, output_stream.str());
}

} // namespace
