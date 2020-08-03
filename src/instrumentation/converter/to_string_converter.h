#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_TO_STRING_CONVERTER_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_TO_STRING_CONVERTER_H_

#include "instrumentation/converter/converter_i.h"

namespace bart {

namespace instrumentation {

namespace converter {

template <typename InputType>
class ToStringConverter : public ConverterI<InputType, std::string> {
 public:
  virtual ~ToStringConverter() = default;
};

} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_TO_STRING_CONVERTER_H_
