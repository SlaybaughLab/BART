#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_FACTORY_HPP_
#define BART_SRC_INSTRUMENTATION_CONVERTER_FACTORY_HPP_

#include "utility/factory/auto_registering_factory.hpp"

namespace bart::instrumentation::converter {

template <typename InputType, typename OutputType> class ConverterI;

enum class ConverterName {
  kCalculatorVectorSubtractor = 0, // calculator/vector_subtractor.h
  kFourierTransform = 1, // fourier/fourier_transform.h
  kConvergenceToString = 2, // to_string/convergence_to_string.h
  kDoubleToString = 3, // to_string/double_to_string.h
  kIntDoublePairToString = 4, // to_string/int_double_pair_to_string.h
  kIntVectorComplexPairToString = 5, // to_string/int_vector_complex_pair_to_string.h
  kStringColorPairToString = 6, // to_string/string_color_pair_to_string.h
  kPairIncrementer = 7,
  kDealiiToComplexVector = 8,
  kGroupScalarFluxExtractor = 9
};

template <typename InputType, typename OutputType, typename ...T>
class ConverterIFactory : public utility::factory::AutoRegisteringFactory<
    ConverterName,
    std::unique_ptr<ConverterI<InputType, OutputType>>(*)(T...)> {};

[[nodiscard]] inline auto to_string(ConverterName to_convert) -> std::string {
  switch (to_convert) {
    case (ConverterName::kCalculatorVectorSubtractor):
      return std::string{"ConverterName::kCalculatorVectorSubtractor"};
    case (ConverterName::kFourierTransform):
      return std::string{"ConverterName::kFourierTransform"};
    case (ConverterName::kConvergenceToString):
      return std::string{"ConverterName::kConvergenceToString"};
    case (ConverterName::kDoubleToString):
      return std::string{"ConverterName::kDoubleToString"};
    case (ConverterName::kIntDoublePairToString):
      return std::string{"ConverterName::kIntDoublePairToString"};
    case (ConverterName::kIntVectorComplexPairToString):
      return std::string{"ConverterName::kIntVectorComplexPairToString"};
    case (ConverterName::kStringColorPairToString):
      return std::string{"ConverterName::kStringColorPairToString"};
    case (ConverterName::kPairIncrementer):
      return std::string{"ConverterName::kPairIncrementer"};
    case (ConverterName::kDealiiToComplexVector):
      return std::string{"ConverterName::kDealiiToComplexVector"};
    case (ConverterName::kGroupScalarFluxExtractor):
      return std::string{"ConverterName::kGroupScalarFluxExtractor"};
  }
  return std::string{"Unknown ConverterName conversion to string requested"};
}

} // namespace bart::instrumentation::converter

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_FACTORY_HPP_
