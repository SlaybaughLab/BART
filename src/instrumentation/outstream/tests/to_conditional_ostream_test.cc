#include "instrumentation/outstream/to_conditional_ostream.h"

#include <deal.II/base/conditional_ostream.h>
#include <memory>
#include <sstream>

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class InstrumentationOutstreamToConditionalOstreamTest : public ::testing::Test {
 public:
  using OutstreamType = instrumentation::outstream::ToConditionalOstream;
  using ConditionalOstreamType = dealii::ConditionalOStream;
  std::unique_ptr<OutstreamType> test_outstream;
  std::ostringstream output_stream;
  void SetUp() override;
};

void InstrumentationOutstreamToConditionalOstreamTest::SetUp() {
  test_outstream = std::make_unique<OutstreamType>(
      std::move(std::make_unique<ConditionalOstreamType>(output_stream, true)));
}

TEST_F(InstrumentationOutstreamToConditionalOstreamTest, Constructor) {
  std::ostringstream string_stream;
  EXPECT_NO_THROW({
    instrumentation::outstream::ToConditionalOstream test(std::move(
        std::make_unique<dealii::ConditionalOStream>(string_stream, true)));
    });
}

TEST_F(InstrumentationOutstreamToConditionalOstreamTest, StreamGetter) {
  EXPECT_NE(nullptr, test_outstream->conditional_ostream_ptr());
}

TEST_F(InstrumentationOutstreamToConditionalOstreamTest, OutstreamString) {
  std::string test_string{"this is a test string\n"},
      second_string{"so is this"};
  test_outstream->Output(test_string).Output(second_string);
  EXPECT_EQ(test_string + second_string, output_stream.str());
}

} // namespace
