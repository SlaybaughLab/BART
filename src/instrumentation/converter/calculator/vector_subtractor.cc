#include "instrumentation/converter/calculator/vector_subtractor.h"
#include "instrumentation/converter/factory.h"

namespace bart {

namespace instrumentation {

namespace converter {

namespace calculator {

bool VectorSubtractor::is_registered_ = ConverterIFactory<DealiiVector, DealiiVector, DealiiVector, AbsoluteValue>::get()
    .RegisterConstructor(
        ConverterName::kCalculatorVectorSubtractor,
        [](DealiiVector minuend, AbsoluteValue calculate_absolute_value) {
          std::unique_ptr<ConverterI<DealiiVector, DealiiVector>> return_ptr =
              std::make_unique<VectorSubtractor>(minuend, calculate_absolute_value);
          return return_ptr;
        });

} // namespace calculator

} // namespace converter

} // namespace instrumentation

} // namespace bart
