#include "instrumentation/instrument.h"

#include "instrumentation/converter/tests/converter_mock.h"
#include "instrumentation/output/tests/output_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class InstrumentationInstrumentTest : public ::testing::Test {
 public:
};

TEST_F(InstrumentationInstrumentTest, Dummy) {
  EXPECT_TRUE(false);
}

} // namespace
