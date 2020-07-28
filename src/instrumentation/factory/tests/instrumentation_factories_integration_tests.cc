#include "instrumentation/factory/instrumentation_factories.h"

#include <deal.II/base/conditional_ostream.h>

#include "instrumentation/outstream/to_conditional_ostream.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

namespace instrumentation = bart::instrumentation;

class InstrumentationFactoriesIntegrationTests : public ::testing::Test {
 public:
};

TEST_F(InstrumentationFactoriesIntegrationTests, ConditionalOstream) {
  auto outstream_ptr = instrumentation::factory::MakeOutstream<std::string>(
      std::make_unique<dealii::ConditionalOStream>(std::cout, false));
  ASSERT_NE(outstream_ptr, nullptr);
  using ExpectedType = instrumentation::outstream::ToConditionalOstream;
  ASSERT_NE(dynamic_cast<ExpectedType*>(outstream_ptr.get()), nullptr);
}

} // namespace
