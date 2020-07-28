#include "instrumentation/factory/instrumentation_factories.h"

#include <deal.II/base/conditional_ostream.h>

#include "instrumentation/output/to_conditional_ostream.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

namespace instrumentation = bart::instrumentation;

class InstrumentationFactoriesIntegrationTests : public ::testing::Test {
 public:
};

TEST_F(InstrumentationFactoriesIntegrationTests, ConditionalOstream) {
  auto outputter_ptr = instrumentation::factory::MakeOutputter<std::string>(
      std::make_unique<dealii::ConditionalOStream>(std::cout, false));
  ASSERT_NE(outputter_ptr, nullptr);
  using ExpectedType = instrumentation::output::ToConditionalOstream;
  ASSERT_NE(dynamic_cast<ExpectedType*>(outputter_ptr.get()), nullptr);
}

} // namespace
