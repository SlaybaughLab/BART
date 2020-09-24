#include "calculator/fourier/fourier_transform_fftw.h"
#include "instrumentation/factory/instrument_factories.h"
#include "instrumentation/factory/component_factories.h"
#include "instrumentation/converter/pair_incrementer.h"
#include "instrumentation/converter/calculator/vector_subtractor.h"
#include "instrumentation/converter/dealii_to_complex_vector.h"
#include "instrumentation/converter/fourier/fourier_transform.h"
#include "instrumentation/converter/multi_converter.h"
#include "instrumentation/converter/to_string/int_vector_complex_pair_to_string.h"
#include "instrumentation/outstream/to_ostream.h"

namespace bart {

namespace instrumentation {

namespace factory {

std::unique_ptr<FourierInstrumentType> MakeErrorFourierTransformInstrument(
    DealiiVector solution, std::unique_ptr<std::ostream> output_stream) {
  // The first converter is a vector subtractor to subtract the input from the
  // solution
  using AbsoluteValue = converter::calculator::AbsoluteValue;
  using ComplexVector = std::vector<std::complex<double>>;
  using VectorSubtractor = converter::calculator::VectorSubtractor;
  using DealiiVectorToComplexVector = converter::DealiiToComplexVector;
  // The first multi-converter takes an input DealiiVector, subtracts it from
  // our solution and converts to a complex vector
  auto complex_error_converter = std::make_unique<
      converter::MultiConverter<DealiiVector, DealiiVector, ComplexVector>>(
          std::make_unique<VectorSubtractor>(solution, AbsoluteValue(false)),
          std::make_unique<DealiiVectorToComplexVector>());
  // The next multi-converter uses that multi-converter, and calculates the
  // fourier transform.
  using FourierTransformCalculator = bart::calculator::fourier::FourierTransformFFTW;
  using FourierTransformConverter = converter::fourier::FourierTransform;
  auto error_fourier_transform_converter = std::make_unique<
      converter::MultiConverter<DealiiVector, ComplexVector, ComplexVector>>(
          std::move(complex_error_converter),
          std::make_unique<FourierTransformConverter>(
              std::make_unique<FourierTransformCalculator>(solution.size())));
  // The next multi-converter changes the complex vector into a pair with an
  // incrementing value
  using PairIncrementer = converter::PairIncrementer<ComplexVector>;
  auto incremented_error_fourier_transform_converter = std::make_unique<
      converter::MultiConverter<DealiiVector, ComplexVector, std::pair<int, ComplexVector>>>(
          std::move(error_fourier_transform_converter),
          std::make_unique<PairIncrementer>());
  // The final multi-converter changes the pair into a string to be read out
  using IntVectorComplexPairToString = converter::to_string::IntVectorComplexPairToString;
  auto final_converter = std::make_unique<
      converter::MultiConverter<DealiiVector, std::pair<int, ComplexVector>, std::string>>(
          std::move(incremented_error_fourier_transform_converter),
          std::make_unique<IntVectorComplexPairToString>());
  // Finally, we build the instrument
  using ReturnType = instrumentation::Instrument<DealiiVector, std::string>;
  return std::make_unique<ReturnType>(
      std::move(final_converter),
      std::make_unique<outstream::ToOstream>(std::move(output_stream)));
}

} // namespace factory

} // namespace instrumentation

} // namespace bart

