#include "instrumentation/factory/instrumentation_factories.h"

#include <deal.II/base/conditional_ostream.h>

#include "instrumentation/converter/convergence_to_string.h"
#include "instrumentation/outstream/to_conditional_ostream.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

namespace instrumentation = bart::instrumentation;
namespace convergence = bart::convergence;

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

TEST_F(InstrumentationFactoriesIntegrationTests, ConverterConvergenceToString) {
  auto converter_ptr =  instrumentation::factory::MakeConverter<convergence::Status, std::string>();
  ASSERT_NE(converter_ptr, nullptr);
  using ExpectedType = instrumentation::converter::ConvergenceToString;
  ASSERT_NE(dynamic_cast<ExpectedType*>(converter_ptr.get()), nullptr);
}

} // namespace
