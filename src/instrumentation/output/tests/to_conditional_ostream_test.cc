#include "instrumentation/output/to_conditional_ostream.h"

#include <deal.II/base/conditional_ostream.h>

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class InstrumentationOutputToConditionalOstreamTest : public ::testing::Test {
 public:
};

TEST_F(InstrumentationOutputToConditionalOstreamTest, Constructor) {
  EXPECT_NO_THROW({
    instrumentation::output::ToConditionalOstream<std::string> test;
  });
}

} // namespace
