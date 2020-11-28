#include "framework/builder/framework_builder.hpp"

#include "convergence/status.hpp"
#include "instrumentation/tests/instrument_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "utility/colors.hpp"

#include <string>

namespace  {

using namespace bart;

class FrameworkBuilderInstrumentationTest : public ::testing::Test {
 public:
  static constexpr int dim { 2 };
  using ColorStatusPair = std::pair<std::string, utility::Color>;
  using ColorStatusInstrument = instrumentation::InstrumentMock<ColorStatusPair>;
  using ConvergenceInstrument = instrumentation::InstrumentMock<convergence::Status>;
  using StatusInstrument = instrumentation::InstrumentMock<std::string>;

  std::shared_ptr<ColorStatusInstrument> color_status_instrument_ptr_{ nullptr };
  std::shared_ptr<ConvergenceInstrument> convergence_instrument_ptr_{ nullptr };
  std::shared_ptr<StatusInstrument> status_instrument_ptr_{ nullptr };

  auto SetUp() -> void override;
};

auto FrameworkBuilderInstrumentationTest::SetUp() -> void {
  color_status_instrument_ptr_ = std::make_shared<ColorStatusInstrument>();
  convergence_instrument_ptr_ = std::make_shared<ConvergenceInstrument>();
  status_instrument_ptr_ = std::make_shared<StatusInstrument>();
}

TEST_F(FrameworkBuilderInstrumentationTest, GettersAndSetters) {
  framework::builder::FrameworkBuilder<dim> test_builder;
  EXPECT_NO_THROW({
    test_builder.set_color_status_instrument_ptr(color_status_instrument_ptr_);
  });
  EXPECT_NO_THROW({
    test_builder.set_convergence_status_instrument_ptr(convergence_instrument_ptr_);
  });
  EXPECT_NO_THROW({
    test_builder.set_status_instrument_ptr(status_instrument_ptr_);
  });

  EXPECT_EQ(test_builder.color_status_instrument_ptr(), color_status_instrument_ptr_);
  EXPECT_EQ(test_builder.status_instrument_ptr(), status_instrument_ptr_);
  EXPECT_EQ(test_builder.convergence_status_instrument_ptr(), convergence_instrument_ptr_);
}

TEST_F(FrameworkBuilderInstrumentationTest, GettersAndSettersNullptrs) {
  framework::builder::FrameworkBuilder<dim> test_builder;
  EXPECT_ANY_THROW({
                    test_builder.set_color_status_instrument_ptr(nullptr);
                  });
  EXPECT_ANY_THROW({
                    test_builder.set_convergence_status_instrument_ptr(nullptr);
                  });
  EXPECT_ANY_THROW({
                    test_builder.set_status_instrument_ptr(nullptr);
                  });
}

} // namespace