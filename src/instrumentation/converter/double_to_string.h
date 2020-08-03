#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_DOUBLE_TO_STRING_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_DOUBLE_TO_STRING_H_

#include <map>
#include <string>

#include "instrumentation/converter/converter_i.h"

namespace bart {

namespace instrumentation {

namespace converter {

class DoubleToString : public ConverterI<double, std::string> {
 public:
  enum OutputTerm {
    kValue,
  };
  std::string Convert(const double &input) const override {
    return std::string();
  }
  virtual ~DoubleToString() = default;

};

} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_DOUBLE_TO_STRING_H_
