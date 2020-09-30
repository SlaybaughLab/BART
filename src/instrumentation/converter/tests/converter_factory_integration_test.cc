#include "instrumentation/converter/factory.h"

#include "instrumentation/converter/calculator/vector_subtractor.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

namespace converter = bart::instrumentation::converter;

TEST(InstrumentationConverterIFactoryTest, VectorSubtractorInstantiation) {
  using ExpectedType = converter::calculator::VectorSubtractor;
  using DealiiVector = dealii::Vector<double>;
  using AbsoluteValue = converter::calculator::AbsoluteValue;

  DealiiVector test_vector(10);

  auto vector_subtractor_ptr = converter::ConverterIFactory<DealiiVector, DealiiVector, DealiiVector, AbsoluteValue>::get()
      .GetConstructor(converter::ConverterName::kCalculatorVectorSubtractor)
          (test_vector, AbsoluteValue(true));
  ASSERT_NE(vector_subtractor_ptr, nullptr);
  ASSERT_NE(dynamic_cast<ExpectedType*>(vector_subtractor_ptr.get()), nullptr);
}

} // namespace
