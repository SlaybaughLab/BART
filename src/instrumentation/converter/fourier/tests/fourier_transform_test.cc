#include "instrumentation/converter/fourier/fourier_transform.h"

#include "calculator/fourier/tests/fourier_transform_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

using namespace bart;

using ::testing::TypedEq;
using ::testing::ContainerEq;
using ::testing::WhenDynamicCastTo, ::testing::NotNull, ::testing::Ref;
using ::testing::Return;

class InstrumentationConverterFourierTransform : public ::testing::Test {
 public:
  using FourierCalculator = calculator::fourier::FourierTransformMock;
  using Normalized = calculator::fourier::Normalized;
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

TEST_F(InstrumentationConverterFourierTransform, Convert) {
  const int vector_size{test_helpers::RandomInt(3, 5)};
  using ComplexVector = std::vector<std::complex<double>>;
  ComplexVector input_vector, output_vector;

  auto random_complex = []{
    return std::complex<double>(test_helpers::RandomDouble(-100, 100),
                                test_helpers::RandomDouble(-100, 100));
  };

  for (int i = 0; i < vector_size; ++i) {
    input_vector.push_back(random_complex());
    output_vector.push_back(random_complex());
  }

  EXPECT_CALL(*fourier_calculator_obs_ptr_,
              CalculateDFT(
                  ::testing::TypedEq<const ComplexVector&>(input_vector),
                  Normalized(true)))
      .WillOnce(Return(output_vector));

  auto conversion_result = test_fourier_converter_->Convert(input_vector);
  EXPECT_THAT(conversion_result, ContainerEq(output_vector));
}

} // namespace
