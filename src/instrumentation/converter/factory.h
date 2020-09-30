#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_FACTORY_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_FACTORY_H_

#include "utility/factory/auto_registering_factory.h"

namespace bart {

namespace instrumentation {

namespace converter {

template <typename InputType, typename OutputType> class ConverterI;

enum class ConverterName {
  kCalculatorVectorSubtractor, // calculator/vector_subtractor.h
  kFourierTransform, // fourier/fourier_transform.h
  kConvergenceToString, // to_string/convergence_to_string.h
  kDoubleToString, // to_string/double_to_string.h
  kIntDoublePairToString, // to_string/int_double_pair_to_string.h
  kIntVectorComplexPairToString, // to_string/int_vector_complex_pair_to_string.h
  kStringColorPairToString, // to_string/string_color_pair_to_string.h
};

template <typename InputType, typename OutputType, typename ...T>
class ConverterIFactory : public utility::factory::AutoRegisteringFactory<
    ConverterName,
    std::unique_ptr<ConverterI<InputType, OutputType>>(*)(T...)> {};

} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_FACTORY_H_
