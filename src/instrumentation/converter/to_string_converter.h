#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_TO_STRING_CONVERTER_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_TO_STRING_CONVERTER_H_

#include "instrumentation/converter/converter_i.h"

#include <map>
#include <string>

namespace bart {

namespace instrumentation {

namespace converter {

template <typename InputType, typename OutputTerm>
class ToStringConverter : public ConverterI<InputType, std::string> {
 public:
  using OutputTermToStringMap = std::map<OutputTerm, std::string>;
  ToStringConverter(
      std::string output_format,
      ToStringConverter::OutputTermToStringMap output_term_to_string_map)
      : output_format_(output_format),
        output_term_to_string_map_(output_term_to_string_map){}
  virtual ~ToStringConverter() = default;

  std::string output_format() const { return output_format_; }
  OutputTermToStringMap output_term_to_string_map() const {
    return output_term_to_string_map_;
  }
 protected:
  std::string output_format_;
  OutputTermToStringMap output_term_to_string_map_;
};

} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_TO_STRING_CONVERTER_H_
