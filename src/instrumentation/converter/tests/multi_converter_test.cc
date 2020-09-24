#include "instrumentation/converter/multi_converter.h"

#include "instrumentation/converter/tests/converter_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

namespace instrumentation = bart::instrumentation;

class InstrumentationConverterMultiConverterTest : public ::testing::Test {
 public:
};

TEST_F(InstrumentationConverterMultiConverterTest, Dummy) {
  EXPECT_TRUE(false);
}

} // namespace