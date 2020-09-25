#include "instrumentation/factory/converter_factory.h"
#include "instrumentation/converter/to_string/double_to_string.h"
#include "instrumentation/factory/converter_types.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

namespace converter = bart::instrumentation::converter;
namespace factory = bart::instrumentation::factory;

class InstrumentationFactoryConverterFactoryTest : public ::testing::Test {
 public:
};

TEST_F(InstrumentationFactoryConverterFactoryTest, Dummy) {
  using ExpectedType = converter::to_string::DoubleToString;
  auto converter_ptr = factory::ConverterFactory<double, std::string>::get().get_generator(
      factory::ConverterName::kDoubleToString)();
  ASSERT_NE(converter_ptr, nullptr);
  EXPECT_NE(dynamic_cast<ExpectedType*>(converter_ptr.get()), nullptr);
}

} // namespace
