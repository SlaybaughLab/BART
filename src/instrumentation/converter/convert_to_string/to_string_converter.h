#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_TO_STRING_CONVERTER_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_TO_STRING_CONVERTER_H_

#include "instrumentation/converter/converter_i.h"

#include <map>
#include <string>
#include <sstream>
#include <vector>
#include <variant>

namespace bart {

namespace instrumentation {

namespace converter {

namespace convert_to_string {

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

  std::string SetOutputFormat(
      const std::vector<std::variant<OutputTerm, std::string>> output_format_parts) {
    std::ostringstream new_format_stream;
    for (const auto output_part : output_format_parts) {
      try {
        new_format_stream << output_term_to_string_map_.at(
            std::get<OutputTerm>(output_part));
      } catch (const std::bad_variant_access&) {
        new_format_stream << std::get<std::string>(output_part);
      }
    }
    output_format_ = new_format_stream.str();
    return output_format_;
  };
  std::string output_format() const { return output_format_; }
  OutputTermToStringMap output_term_to_string_map() const {
    return output_term_to_string_map_;
  }
 protected:
  std::string output_format_;
  OutputTermToStringMap output_term_to_string_map_;
};

} // namespace convert_to_string

} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_TO_STRING_CONVERTER_H_
