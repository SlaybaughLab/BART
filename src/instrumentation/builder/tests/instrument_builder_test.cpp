#include "instrumentation/builder/instrument_builder.hpp"

#include <string>
#include <utility>

#include <deal.II/lac/vector.h>

#include "convergence/status.h"
#include "instrumentation/converter/dealii_to_complex_vector.h"
#include "instrumentation/converter/calculator/vector_subtractor.h"
#include "instrumentation/converter/to_string/string_color_pair_to_string.h"
#include "instrumentation/converter/to_string/convergence_to_string.h"
#include "instrumentation/converter/to_string/int_double_pair_to_string.h"
#include "instrumentation/converter/to_string/int_vector_complex_pair_to_string.h"
#include "instrumentation/converter/multi_converter.hpp"
#include "instrumentation/converter/pair_incrementer.h"
#include "instrumentation/converter/fourier/fourier_transform.h"
#include "instrumentation/converter/system/group_scalar_flux_extractor.hpp"
#include "instrumentation/outstream/to_conditional_ostream.h"
#include "instrumentation/outstream/to_ostream.h"
#include "instrumentation/instrument.h"
#include "instrumentation/basic_instrument.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"
#include "system/moments/spherical_harmonic_i.h"
#include "utility/colors.h"

namespace  {

namespace instrumentation = bart::instrumentation;
namespace test_helpers = bart::test_helpers;

class InstrumentationBuilderInstrumentBuilderTest : public ::testing::Test {
 public:
  using Builder = instrumentation::builder::InstrumentBuilder;
  using DealiiVector = dealii::Vector<double>;
  using InstrumentName = instrumentation::builder::InstrumentName;
  using IntDoublePair = std::pair<int, double>;
  using ConvergenceStatus = bart::convergence::Status;
  using StringColorPair = std::pair<std::string, bart::utility::Color>;
  using SphericalHarmonics = bart::system::moments::SphericalHarmonicI;
};

TEST_F(InstrumentationBuilderInstrumentBuilderTest,
       ColorStringToConditionalOstream) {
  using InstrumentType = instrumentation::Instrument<StringColorPair, std::string>;
  using ConverterType = instrumentation::converter::to_string::StringColorPairToString;
  using OutstreamType = instrumentation::outstream::ToConditionalOstream;
  auto instrument_ptr = Builder::BuildInstrument<StringColorPair>(
      instrumentation::builder::InstrumentName::kColorStatusToConditionalOstream);
  ASSERT_NE(instrument_ptr, nullptr);
  auto dynamic_ptr = dynamic_cast<InstrumentType*>(instrument_ptr.get());
  ASSERT_NE(dynamic_ptr, nullptr);
  ASSERT_NE(dynamic_cast<ConverterType*>(dynamic_ptr->converter_ptr()), nullptr);
  ASSERT_NE(dynamic_cast<OutstreamType*>(dynamic_ptr->outstream_ptr()), nullptr);
}

TEST_F(InstrumentationBuilderInstrumentBuilderTest, FourierInstrument) {
  using InstrumentType = instrumentation::Instrument<SphericalHarmonics, std::string>;
  using OutStreamType = instrumentation::outstream::ToOstream;
  const std::string filename{ "filename.csv" };
  const int group{ test_helpers::RandomInt(0, 10) };

  dealii::Vector<double> error_vector(test_helpers::RandomInt(5, 10));
  for (auto& entry : error_vector)
    entry = test_helpers::RandomDouble(0, 1000);

  auto instrument_ptr = Builder::BuildInstrument<SphericalHarmonics>(
      InstrumentName::kFourierOfScalarFluxErrorToFile, group, error_vector, filename);
  ASSERT_NE(instrument_ptr, nullptr);
  auto dynamic_ptr = dynamic_cast<InstrumentType*>(instrument_ptr.get());
  ASSERT_NE(dynamic_ptr, nullptr);

  // This is a nested instrument with many levels; we'll check each converter
  using ComplexVector = std::vector<std::complex<double>>;
  using IntComplexVectorPair = std::pair<int, ComplexVector>;
  using MultiConverterOne = instrumentation::converter::MultiConverter<SphericalHarmonics, IntComplexVectorPair, std::string>;
  // First multi-converter
  auto multi_converter_one_ptr =
      dynamic_cast<MultiConverterOne*>(dynamic_ptr->converter_ptr());
  ASSERT_NE(multi_converter_one_ptr, nullptr);
  ASSERT_NE(
      dynamic_cast<instrumentation::converter::to_string::IntVectorComplexPairToString*>(
          multi_converter_one_ptr->second_stage_converter_ptr()), nullptr);
  // Next multi-converter
  using MultiConverterTwo = instrumentation::converter::MultiConverter<SphericalHarmonics, ComplexVector, IntComplexVectorPair>;
  auto multi_converter_two_ptr = dynamic_cast<MultiConverterTwo*>(
      multi_converter_one_ptr->first_stage_converter_ptr());
  ASSERT_NE(multi_converter_two_ptr, nullptr);
  ASSERT_NE(
      dynamic_cast<instrumentation::converter::PairIncrementer<ComplexVector>*>(
          multi_converter_two_ptr->second_stage_converter_ptr()), nullptr);
  // Next multi-converter
  using MultiConverterThree = instrumentation::converter::MultiConverter<SphericalHarmonics, ComplexVector, ComplexVector>;
  auto multi_converter_three_ptr = dynamic_cast<MultiConverterThree*>(
      multi_converter_two_ptr->first_stage_converter_ptr());
  ASSERT_NE(multi_converter_three_ptr, nullptr);
  ASSERT_NE(
      dynamic_cast<instrumentation::converter::fourier::FourierTransform*>(
          multi_converter_three_ptr->second_stage_converter_ptr()), nullptr);
  // Next multi-converter
  using MultiConverterFour = instrumentation::converter::MultiConverter<SphericalHarmonics , DealiiVector , ComplexVector >;
  auto multi_converter_four_ptr = dynamic_cast<MultiConverterFour*>(
      multi_converter_three_ptr->first_stage_converter_ptr());
  ASSERT_NE(multi_converter_four_ptr, nullptr);
  ASSERT_NE(dynamic_cast<instrumentation::converter::DealiiToComplexVector*>(
      multi_converter_four_ptr->second_stage_converter_ptr()), nullptr);
  ASSERT_NE(dynamic_cast<OutStreamType*>(dynamic_ptr->outstream_ptr()), nullptr);
  EXPECT_EQ(remove(filename.c_str()), 0);

  // Next multi-converter
  using MultiConverterFive = instrumentation::converter::MultiConverter<SphericalHarmonics , DealiiVector , DealiiVector >;
  auto multi_converter_five_ptr = dynamic_cast<MultiConverterFive*>(
      multi_converter_four_ptr->first_stage_converter_ptr());
  ASSERT_NE(multi_converter_five_ptr, nullptr);

  auto group_scalar_flux_extractor_ptr =
      dynamic_cast<instrumentation::converter::system::GroupScalarFluxExtractor*>(
          multi_converter_five_ptr->first_stage_converter_ptr());
  ASSERT_NE(group_scalar_flux_extractor_ptr, nullptr);
  EXPECT_EQ(group_scalar_flux_extractor_ptr->group_to_extract(), group);

  auto vector_subtractor_ptr = dynamic_cast<instrumentation::converter::calculator::VectorSubtractor*>(
      multi_converter_five_ptr->second_stage_converter_ptr());
  ASSERT_NE(vector_subtractor_ptr, nullptr);
  EXPECT_EQ(vector_subtractor_ptr->minuend(), error_vector);

}

TEST_F(InstrumentationBuilderInstrumentBuilderTest,
       ConvergenceStatusToConditionalOstream) {
  using InstrumentType = instrumentation::Instrument<ConvergenceStatus, std::string>;
  using ConverterType = instrumentation::converter::to_string::ConvergenceToString;
  using OutstreamType = instrumentation::outstream::ToConditionalOstream;
  auto instrument_ptr = Builder::BuildInstrument<ConvergenceStatus>(
      InstrumentName::kConvergenceStatusToConditionalOstream);
  ASSERT_NE(instrument_ptr, nullptr);
  auto dynamic_ptr = dynamic_cast<InstrumentType*>(instrument_ptr.get());
  ASSERT_NE(dynamic_ptr, nullptr);
  ASSERT_NE(dynamic_cast<ConverterType*>(dynamic_ptr->converter_ptr()), nullptr);
  ASSERT_NE(dynamic_cast<OutstreamType*>(dynamic_ptr->outstream_ptr()), nullptr);
}

TEST_F(InstrumentationBuilderInstrumentBuilderTest,
       IntDoublePairToFile) {
  const std::string filename{ "filename.csv" };
  using InstrumentType = instrumentation::Instrument<IntDoublePair, std::string>;
  using ConverterType = instrumentation::converter::to_string::IntDoublePairToString;
  using OutstreamType = instrumentation::outstream::ToOstream;
  auto instrument_ptr = Builder::BuildInstrument<IntDoublePair>(
      InstrumentName::kIntDoublePairToFile,
      filename);
  ASSERT_NE(instrument_ptr, nullptr);
  auto dynamic_ptr = dynamic_cast<InstrumentType*>(instrument_ptr.get());
  ASSERT_NE(dynamic_ptr, nullptr);
  ASSERT_NE(dynamic_cast<ConverterType*>(dynamic_ptr->converter_ptr()), nullptr);
  ASSERT_NE(dynamic_cast<OutstreamType*>(dynamic_ptr->outstream_ptr()), nullptr);
  EXPECT_EQ(remove(filename.c_str()), 0);
}

TEST_F(InstrumentationBuilderInstrumentBuilderTest,
       StringToConditionalOstream) {
  using InstrumentType = instrumentation::BasicInstrument<std::string>;
  using OutstreamType = instrumentation::outstream::ToConditionalOstream;
  auto instrument_ptr = Builder::BuildInstrument<std::string>(
      instrumentation::builder::InstrumentName::kStringToConditionalOstream);
  ASSERT_NE(instrument_ptr, nullptr);
  auto dynamic_ptr = dynamic_cast<InstrumentType*>(instrument_ptr.get());
  ASSERT_NE(dynamic_ptr, nullptr);
  ASSERT_NE(dynamic_cast<OutstreamType*>(dynamic_ptr->outstream_ptr()), nullptr);
}

TEST_F(InstrumentationBuilderInstrumentBuilderTest,
       BadInstrumentNamesPerType) {
  EXPECT_ANY_THROW(Builder::BuildInstrument<std::string>(
      instrumentation::builder::InstrumentName::kColorStatusToConditionalOstream));
  EXPECT_ANY_THROW(Builder::BuildInstrument<StringColorPair>(
      instrumentation::builder::InstrumentName::kStringToConditionalOstream));

  for (auto bad_name : {InstrumentName::kColorStatusToConditionalOstream,
                        InstrumentName::kStringToConditionalOstream}) {
    EXPECT_ANY_THROW(Builder::BuildInstrument<ConvergenceStatus>(bad_name));
  }

  for (auto bad_name : {InstrumentName::kColorStatusToConditionalOstream,
                        InstrumentName::kStringToConditionalOstream,
                        InstrumentName::kConvergenceStatusToConditionalOstream}) {
    EXPECT_ANY_THROW(Builder::BuildInstrument<IntDoublePair>(bad_name, std::string{}));
  }

  for (auto bad_name : {InstrumentName::kColorStatusToConditionalOstream,
                        InstrumentName::kStringToConditionalOstream,
                        InstrumentName::kConvergenceStatusToConditionalOstream,
                        InstrumentName::kIntDoublePairToFile}) {
    EXPECT_ANY_THROW(Builder::BuildInstrument<SphericalHarmonics>(bad_name,
                                                            int{},
                                                            DealiiVector{},
                                                            std::string{}));
  }
}



} // namespace
