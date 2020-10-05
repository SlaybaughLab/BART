#include "instrumentation/outstream/factory.h"

#include <memory>
#include <sstream>

#include "instrumentation/outstream/to_ostream.h"
#include "instrumentation/outstream/to_conditional_ostream.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

namespace outstream = bart::instrumentation::outstream;

class InstrumentationOutstreamIFactoryTest : public ::testing::Test {
 public:
  using OutstreamName = outstream::OutstreamName;
};

TEST_F(InstrumentationOutstreamIFactoryTest, ToOstreamInstantiation) {
  using ExpectedType = outstream::ToOstream;
  std::unique_ptr<std::ostringstream> out_string_stream;
  auto outstream_ptr =
      outstream::OutstreamIFactory<std::string, std::unique_ptr<std::ostream>>::get()
      .GetConstructor(OutstreamName::kToOstream)(std::move(out_string_stream));
  ASSERT_NE(outstream_ptr, nullptr);
  ASSERT_NE(dynamic_cast<ExpectedType*>(outstream_ptr.get()), nullptr);
}

TEST_F(InstrumentationOutstreamIFactoryTest, ToConditionalOstream) {
  using ExpectedType = outstream::ToConditionalOstream;
  using ConditionalOstreamPtrType = std::unique_ptr<dealii::ConditionalOStream>;
  auto pout_ptr = std::make_unique<dealii::ConditionalOStream>(std::cout);
  auto outstream_ptr =
      outstream::OutstreamIFactory<std::string, ConditionalOstreamPtrType>::get()
      .GetConstructor(OutstreamName::kToConditionalOstream)(std::move(pout_ptr));
  ASSERT_NE(outstream_ptr, nullptr);
  ASSERT_NE(dynamic_cast<ExpectedType*>(outstream_ptr.get()), nullptr);
}

} // namespace
