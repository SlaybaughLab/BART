#include "instrumentation/converter/dealii_to_complex_vector.h"

namespace bart::instrumentation::converter {

bool DealiiToComplexVector::is_registered_ =
    ConverterIFactory<DealiiVector, ComplexVector>::get()
        .RegisterConstructor(
            ConverterName::kDealiiToComplexVector,
            []() {
              std::unique_ptr<ConverterI<DealiiVector, ComplexVector>> return_ptr =
                  std::make_unique<DealiiToComplexVector>();
              return return_ptr;
            });

} // namespace bart::instrumentation::converter
