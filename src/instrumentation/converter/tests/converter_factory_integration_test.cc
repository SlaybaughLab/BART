#include "instrumentation/converter/factory.hpp"

#include "convergence/status.h"
#include "calculator/fourier/tests/fourier_transform_mock.h"
#include "instrumentation/converter/calculator/vector_subtractor.h"
#include "instrumentation/converter/fourier/fourier_transform.h"
#include "instrumentation/converter/to_string/convergence_to_string.h"
#include "instrumentation/converter/to_string/double_to_string.h"
#include "instrumentation/converter/to_string/int_double_pair_to_string.h"
#include "instrumentation/converter/to_string/int_vector_complex_pair_to_string.h"
#include "instrumentation/converter/to_string/string_color_pair_to_string.h"
#include "instrumentation/converter/pair_incrementer.h"
#include "instrumentation/converter/dealii_to_complex_vector.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

namespace calculator = bart::calculator;
namespace convergence = bart::convergence;
namespace converter = bart::instrumentation::converter;

class InstrumentationConverterIFactoryTest : public ::testing::Test {
 public:
  using AbsoluteValue = converter::calculator::AbsoluteValue;
  using Color = bart::utility::Color;
  using ComplexVector = std::vector<std::complex<double>>;
  using DealiiVector = dealii::Vector<double>;
  using FourierCalculator = calculator::fourier::FourierTransformI;
  using FourierCalculatorMock = calculator::fourier::FourierTransformMock;
  using IntComplexVectorPair = std::pair<int, ComplexVector>;
  using Normalized = calculator::fourier::Normalized;
};

TEST_F(InstrumentationConverterIFactoryTest, VectorSubtractorInstantiation) {
  using ExpectedType = converter::calculator::VectorSubtractor;

  DealiiVector test_vector(10);
  auto vector_subtractor_ptr = converter::ConverterIFactory<DealiiVector, DealiiVector, DealiiVector, AbsoluteValue>::get()
      .GetConstructor(converter::ConverterName::kCalculatorVectorSubtractor)
          (test_vector, AbsoluteValue(true));
  ASSERT_NE(vector_subtractor_ptr, nullptr);
  ASSERT_NE(dynamic_cast<ExpectedType*>(vector_subtractor_ptr.get()), nullptr);
}

TEST_F(InstrumentationConverterIFactoryTest, FourierTransformInstantiation) {
  using ExpectedType = converter::fourier::FourierTransform;

  auto fourier_transform_ptr = converter::ConverterIFactory<ComplexVector, ComplexVector, std::unique_ptr<FourierCalculator>, Normalized>::get()
      .GetConstructor(converter::ConverterName::kFourierTransform)
          (std::make_unique<FourierCalculatorMock>(), Normalized(true));
  ASSERT_NE(fourier_transform_ptr, nullptr);
  auto dynamic_ptr = dynamic_cast<ExpectedType*>(fourier_transform_ptr.get());
  ASSERT_NE(dynamic_ptr, nullptr);
  EXPECT_TRUE(dynamic_ptr->returns_normalized());
  EXPECT_NE(dynamic_ptr->fourier_calculator_ptr(), nullptr);
}

TEST_F(InstrumentationConverterIFactoryTest, ConvergenceToStringInstantiation) {
  using ExpectedType = converter::to_string::ConvergenceToString;
  auto convergence_to_string_ptr = converter::ConverterIFactory<convergence::Status, std::string>::get()
      .GetConstructor(converter::ConverterName::kConvergenceToString)();
  ASSERT_NE(convergence_to_string_ptr, nullptr);
  ASSERT_NE(dynamic_cast<ExpectedType*>(convergence_to_string_ptr.get()),
            nullptr);
}

TEST_F(InstrumentationConverterIFactoryTest, DoubletoStringInstantiation) {
  using ExpectedType = converter::to_string::DoubleToString;
  auto double_to_string_ptr = converter::ConverterIFactory<double, std::string>::get()
      .GetConstructor(converter::ConverterName::kDoubleToString)();
  ASSERT_NE(double_to_string_ptr, nullptr);
  ASSERT_NE(dynamic_cast<ExpectedType*>(double_to_string_ptr.get()), nullptr);
}

TEST_F(InstrumentationConverterIFactoryTest, IntDoublePairToStringInstantiation) {
  using ExpectedType = converter::to_string::IntDoublePairToString;
  auto converter_ptr = converter::ConverterIFactory<std::pair<int, double>, std::string>::get()
      .GetConstructor(converter::ConverterName::kIntDoublePairToString)();
  ASSERT_NE(converter_ptr, nullptr);
  ASSERT_NE(dynamic_cast<ExpectedType*>(converter_ptr.get()), nullptr);
}

TEST_F(InstrumentationConverterIFactoryTest,
       IntVectorComplexPairToStringInstantiation) {
  using ExpectedType = converter::to_string::IntVectorComplexPairToString;
  auto converter_ptr = converter::ConverterIFactory<IntComplexVectorPair, std::string>::get()
      .GetConstructor(converter::ConverterName::kIntVectorComplexPairToString)();
  ASSERT_NE(converter_ptr, nullptr);
  ASSERT_NE(dynamic_cast<ExpectedType*>(converter_ptr.get()), nullptr);
}

TEST_F(InstrumentationConverterIFactoryTest,
       StringColorPairToStringInstantiation) {
  using ExpectedType = converter::to_string::StringColorPairToString;
  using ColorStringPair = std::pair<std::string, Color>;
  auto converter_ptr = converter::ConverterIFactory<ColorStringPair, std::string>::get()
      .GetConstructor(converter::ConverterName::kStringColorPairToString)();
  ASSERT_NE(converter_ptr, nullptr);
  ASSERT_NE(dynamic_cast<ExpectedType*>(converter_ptr.get()), nullptr);
}

TEST_F(InstrumentationConverterIFactoryTest,
       PairIncrementerInstantiation) {
  using ExpectedType = converter::PairIncrementer<std::string>;

  auto converter_ptr = converter::ConverterIFactory<std::string, std::pair<int, std::string>>::get()
      .GetConstructor(converter::ConverterName::kPairIncrementer)();
  ASSERT_NE(converter_ptr, nullptr);
  ASSERT_NE(dynamic_cast<ExpectedType*>(converter_ptr.get()), nullptr);
}

TEST_F(InstrumentationConverterIFactoryTest,
    DealiiToComplexVectorInstantiation) {
  using ExpectedType = converter::DealiiToComplexVector;

  auto converter_ptr = converter::ConverterIFactory<DealiiVector, ComplexVector>::get()
      .GetConstructor(converter::ConverterName::kDealiiToComplexVector)();
  ASSERT_NE(converter_ptr, nullptr);
  ASSERT_NE(dynamic_cast<ExpectedType*>(converter_ptr.get()), nullptr);
}

} // namespace
