#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_MULTI_CONVERTER_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_MULTI_CONVERTER_H_

#include "instrumentation/converter/converter_i.h"

namespace bart {

namespace instrumentation {

namespace converter {

template <typename InputType, typename IntermediateType, typename OutputType>
class MultiConverter : ConverterI<InputType, OutputType> {
 public:
 private:
  OutputType Convert(const InputType &input) const override {
    return nullptr;
  }
};

} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_MULTI_CONVERTER_H_
