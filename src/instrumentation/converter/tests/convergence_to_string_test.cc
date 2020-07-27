#include "instrumentation/converter/convergence_to_string.h"

#include "convergence/status.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class InstrumentationConverterConvergenceToStringTest : public ::testing::Test {
 public:
  
};

TEST_F(InstrumentationConverterConvergenceToStringTest, Constructor) {
  using ConverterType = instrumentation::converter::ConvergenceToString;
  EXPECT_NO_THROW({
    ConverterType test_converter;
  });
}

} // namespace
