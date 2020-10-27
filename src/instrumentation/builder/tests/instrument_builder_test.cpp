#include "instrumentation/builder/instrument_builder.hpp"

#include <string>
#include <utility>

#include <deal.II/lac/vector.h>

#include "convergence/status.hpp"
#include "instrumentation/converter/dealii_to_complex_vector.h"
#include "instrumentation/converter/calculator/vector_subtractor.h"
#include "instrumentation/converter/convert_to_string/string_color_pair_to_string.h"
#include "instrumentation/converter/convert_to_string/convergence_to_string.h"
#include "instrumentation/converter/convert_to_string/int_double_pair_to_string.h"
#include "instrumentation/converter/convert_to_string/int_vector_complex_pair_to_string.h"
#include "instrumentation/converter/multi_converter.hpp"
#include "instrumentation/converter/pair_incrementer.h"
#include "instrumentation/converter/fourier/fourier_transform.h"
#include "instrumentation/converter/system/group_scalar_flux_extractor.hpp"
#include "instrumentation/outstream/to_conditional_ostream.h"
#include "instrumentation/outstream/to_ostream.hpp"
#include "instrumentation/instrument_array.hpp"
#include "instrumentation/instrument.h"
#include "instrumentation/basic_instrument.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"
#include "system/moments/spherical_harmonic_i.h"
#include "system/moments/tests/spherical_harmonic_mock.h"
#include "utility/colors.h"

namespace  {

namespace instrumentation = bart::instrumentation;
namespace test_helpers = bart::test_helpers;

using ::testing::Return, ::testing::ReturnRef, ::testing::DoDefault;

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
  using ConverterType = instrumentation::converter::convert_to_string::StringColorPairToString;
  using OutstreamType = instrumentation::outstream::ToConditionalOstream;
  auto instrument_ptr = Builder::BuildInstrument<StringColorPair>(
      instrumentation::builder::InstrumentName::kColorStatusToConditionalOstream);
  ASSERT_NE(instrument_ptr, nullptr);
  auto dynamic_ptr = dynamic_cast<InstrumentType*>(instrument_ptr.get());
  ASSERT_NE(dynamic_ptr, nullptr);
  ASSERT_NE(dynamic_cast<ConverterType*>(dynamic_ptr->converter_ptr()), nullptr);
  ASSERT_NE(dynamic_cast<OutstreamType*>(dynamic_ptr->outstream_ptr()), nullptr);
}

TEST_F(InstrumentationBuilderInstrumentBuilderTest, AllGroupsFourierInstrument) {
  using InstrumentType = instrumentation::Instrument<SphericalHarmonics, std::string>;
  using SphericalHarmonicMock = bart::system::moments::SphericalHarmonicMock;
  std::unique_ptr<SphericalHarmonics> mock_spherical_harmonics =
      std::make_unique<SphericalHarmonicMock>();
  SphericalHarmonicMock* mock_spherical_harmonic_obs_ptr =
      dynamic_cast<SphericalHarmonicMock*>(mock_spherical_harmonics.get());
  const int total_groups{ test_helpers::RandomInt(5, 10) };
  const int flux_size { test_helpers::RandomInt(10, 20) };
  const std::string filename_base { "filename_base" };
  std::vector<dealii::Vector<double>> group_scalar_fluxes(total_groups);

  for (int i = 0; i < total_groups; ++i) {
    dealii::Vector<double> group_flux(flux_size);
    for (dealii::Vector<double>::size_type i = 0; i < group_flux.size(); ++i) {
      group_flux[i] = test_helpers::RandomDouble(-100, 100);
    }
    group_scalar_fluxes.emplace_back(group_flux);
  }

  EXPECT_CALL(*mock_spherical_harmonic_obs_ptr, total_groups())
      .WillOnce(Return(total_groups));
  for (int i = 0; i < total_groups; ++i) {
    EXPECT_CALL(*mock_spherical_harmonic_obs_ptr,
                GetMoment(std::array{i, 0, 0}))
        .WillOnce(ReturnRef(group_scalar_fluxes.at(i)));
  }

  auto instrument_ptr = Builder::BuildInstrument<SphericalHarmonics>(
      InstrumentName::kFourierTransformOfAllGroupScalarFluxErrorToFile,
      mock_spherical_harmonics.get(),
      filename_base);

  using InstrumentArrayType = instrumentation::InstrumentArray<SphericalHarmonics>;
  using InstrumentType = instrumentation::Instrument<SphericalHarmonics, std::string>;
  using OutStreamType = instrumentation::outstream::ToOstream;

  auto array_dynamic_ptr = dynamic_cast<InstrumentArrayType*>(instrument_ptr.get());
  ASSERT_NE(array_dynamic_ptr, nullptr);
  EXPECT_EQ(array_dynamic_ptr->size(), total_groups);
  for (const auto& instrument : *array_dynamic_ptr) {
    auto dynamic_ptr = dynamic_cast<InstrumentType*>(instrument.get());
    ASSERT_NE(dynamic_ptr, nullptr);
    ASSERT_NE(dynamic_cast<OutStreamType*>(dynamic_ptr->outstream_ptr()), nullptr);
  }

  for (int i = 0; i < total_groups; ++i) {
    std::string filename = filename_base + '_' + std::to_string(i) + ".csv";
    EXPECT_EQ(remove(filename.c_str()), 0);
  }

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
      InstrumentName::kFourierTransformOfSingleGroupScalarFluxErrorToFile, group, error_vector, filename);
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
      dynamic_cast<instrumentation::converter::convert_to_string::IntVectorComplexPairToString*>(
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
  using ConverterType = instrumentation::converter::convert_to_string::ConvergenceToString;
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
  using ConverterType = instrumentation::converter::convert_to_string::IntDoublePairToString;
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
                        InstrumentName::kIntDoublePairToFile,
                        InstrumentName::kFourierTransformOfAllGroupScalarFluxErrorToFile}) {
    EXPECT_ANY_THROW(Builder::BuildInstrument<SphericalHarmonics>(bad_name,
                                                            int{},
                                                            DealiiVector{},
                                                            std::string{}));
  }
  for (auto bad_name : {InstrumentName::kColorStatusToConditionalOstream,
                        InstrumentName::kStringToConditionalOstream,
                        InstrumentName::kConvergenceStatusToConditionalOstream,
                        InstrumentName::kIntDoublePairToFile,
                        InstrumentName::kFourierTransformOfSingleGroupScalarFluxErrorToFile}) {
    std::unique_ptr<SphericalHarmonics> spherical_harmonics
        = std::make_unique<bart::system::moments::SphericalHarmonicMock>();
    EXPECT_ANY_THROW(Builder::BuildInstrument<SphericalHarmonics>(
        bad_name,
        spherical_harmonics.get(),
        std::string{}));
  }
}



} // namespace
