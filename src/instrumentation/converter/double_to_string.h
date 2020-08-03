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
  enum OutputTerm {kValue, };
  using OutputTermToStringMap = std::map<OutputTerm, std::string>;

  virtual ~DoubleToString() = default;

  std::string Convert(const double &input) const override;

  std::string output_format() const { return output_format_; }
  int precision() const { return precision_; }
  OutputTermToStringMap output_term_to_string_map() const {
    return output_term_to_string_map_; }
 private:
  int precision_ = 2;
  std::string output_format_{"${VALUE}"};
  OutputTermToStringMap output_term_to_string_map_{
      {kValue, "${VALUE}"}
  };
};

} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_DOUBLE_TO_STRING_H_
