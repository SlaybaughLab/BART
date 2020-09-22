#include "instrumentation/converter/fourier/fourier_transform.h"

#include "calculator/fourier/tests/fourier_transform_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

using ::testing::WhenDynamicCastTo, ::testing::NotNull;

class InstrumentationConverterFourierTransform : public ::testing::Test {
 public:
  using FourierCalculator = calculator::fourier::FourierTransformMock;
  using TestConverter = instrumentation::converter::fourier::FourierTransform;
  std::unique_ptr<TestConverter> test_fourier_converter_;

  FourierCalculator* fourier_calculator_obs_ptr_;
  void SetUp() override;
};

void InstrumentationConverterFourierTransform::SetUp() {
  auto fourier_calculator_ptr = std::make_unique<FourierCalculator>();
  fourier_calculator_obs_ptr_ = fourier_calculator_ptr.get();
  test_fourier_converter_ = std::make_unique<TestConverter>(
      std::move(fourier_calculator_ptr));
}

TEST_F(InstrumentationConverterFourierTransform, ConstructorNullDependency) {
  EXPECT_ANY_THROW({
    TestConverter test_converter(nullptr);
  });
}

TEST_F(InstrumentationConverterFourierTransform, DepdendencyGetter) {
  EXPECT_THAT(test_fourier_converter_->fourier_calculator_ptr(),
              WhenDynamicCastTo<FourierCalculator*>(NotNull()));
}

} // namespace
