#include "instrumentation/converter/dealii_to_complex_vector.h"

#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

using namespace bart;
using ::testing::ContainerEq;

class InstrumentationConverterDealiiToComplexVectorTest
    : public ::testing::Test {
 public:
  using TestConverter = instrumentation::converter::DealiiToComplexVector;
  using DealiiVector = instrumentation::converter::DealiiToComplexVector::DealiiVector;
  using ComplexVector = instrumentation::converter::DealiiToComplexVector::ComplexVector;
};

TEST_F(InstrumentationConverterDealiiToComplexVectorTest, Constructor) {
  EXPECT_NO_THROW({
    TestConverter test_converter;
  });
}

TEST_F(InstrumentationConverterDealiiToComplexVectorTest, Convert) {
  const int vector_length{test_helpers::RandomInt(10, 20)};
  DealiiVector input_vector(vector_length);
  ComplexVector expected_output_vector(vector_length);
  for (int i = 0; i < vector_length; ++i) {
    auto entry_value = test_helpers::RandomDouble(-100, 100);
    input_vector[i] = entry_value;
    expected_output_vector.at(i).real(entry_value);
  }
  TestConverter test_converter;
  auto output_vector = test_converter.Convert(input_vector);
  ASSERT_TRUE(output_vector.size() > 0);
  EXPECT_THAT(output_vector, ContainerEq(expected_output_vector));
}

TEST_F(InstrumentationConverterDealiiToComplexVectorTest, ConvertEmpty) {
  TestConverter test_converter;
  auto output_vector = test_converter.Convert(DealiiVector{});
  EXPECT_TRUE(output_vector.size() == 0);
}

} // namespace
