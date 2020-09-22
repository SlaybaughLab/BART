#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_DOUBLE_TO_STRING_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_DOUBLE_TO_STRING_H_

#include <map>
#include <string>
#include <vector>
#include <variant>

#include "instrumentation/converter/to_string/to_string_converter.h"

namespace bart {

namespace instrumentation {

namespace converter {

namespace to_string {

enum class DoubleToStringOutputTerm { kValue = 0, };

class DoubleToString :
    public ToStringConverter<double, DoubleToStringOutputTerm> {
 public:
  using OutputTerm = DoubleToStringOutputTerm;
  using OutputTermToStringMap = std::map<OutputTerm, std::string>;

  DoubleToString();
  virtual ~DoubleToString() = default;

  std::string Convert(const double &input) const override;

  DoubleToString& set_precision(const int to_set) {
    precision_ = to_set;
    return *this; }

  int precision() const { return precision_; }

 private:
  int precision_ = 2;
};

} // namespace to_string

} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_DOUBLE_TO_STRING_H_
