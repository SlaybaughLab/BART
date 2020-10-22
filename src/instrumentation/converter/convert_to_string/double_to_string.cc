#include "instrumentation/converter/convert_to_string/double_to_string.h"

#include "instrumentation/converter/factory.hpp"

#include <iomanip>
#include <sstream>

namespace bart {

namespace instrumentation {

namespace converter {

namespace convert_to_string {

namespace  {
using OutputTerm = DoubleToStringOutputTerm;
std::string default_output_format{"${VALUE}\n"};
std::map<DoubleToStringOutputTerm, std::string>
    default_output_term_to_string_map{{OutputTerm::kValue, "${VALUE}"}};
} // namespace

DoubleToString::DoubleToString()
    : ToStringConverter<double, DoubleToStringOutputTerm>(
    default_output_format, default_output_term_to_string_map) {}

std::string DoubleToString::Convert(const double &input) const {
  std::string return_string{output_format_};
  std::string value_string{output_term_to_string_map_.at(OutputTerm::kValue)};

  if (auto index = return_string.find(value_string);
      index != std::string::npos) {
    std::ostringstream value_stream;
    value_stream << std::fixed << std::setprecision(precision_) << input;
    return_string.replace(index, value_string.size(), value_stream.str());
  }
  return return_string;
}

bool DoubleToString::is_registered_ =
    converter::ConverterIFactory<double, std::string>::get()
    .RegisterConstructor(converter::ConverterName::kDoubleToString,
        [](){
      std::unique_ptr<ConverterI<double, std::string>> return_ptr =
          std::make_unique<DoubleToString>();
          return return_ptr;
        });

} // namespace convert_to_string

} // namespace converter

} // namespace instrumentation

} // namespace bart
