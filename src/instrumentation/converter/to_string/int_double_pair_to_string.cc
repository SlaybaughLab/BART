#include "instrumentation/converter/to_string/int_double_pair_to_string.h"

#include "instrumentation/converter/factory.hpp"

#include <iomanip>

namespace bart {

namespace instrumentation {

namespace converter {

namespace to_string {

namespace  {
using OutputTerm = IntDoublePairToStringOutputTerm;
std::string default_output_format{"${INDEX}, ${VALUE}\n"};
std::map<IntDoublePairToStringOutputTerm, std::string>
    default_output_term_to_string_map{
    {OutputTerm::kValue, "${VALUE}"},
    {OutputTerm::kIndex, "${INDEX}"}};
}

IntDoublePairToString::IntDoublePairToString()
    : ToStringConverter<std::pair<int, double>, OutputTerm>(
    default_output_format, default_output_term_to_string_map){}

std::string IntDoublePairToString::Convert(
    const std::pair<int, double> &input) const {
  std::string return_string{output_format_};
  auto& [index, value] = input;
  std::string index_string{output_term_to_string_map_.at(OutputTerm::kIndex)};
  std::string value_string{output_term_to_string_map_.at(OutputTerm::kValue)};

  if (auto index_index = return_string.find(index_string);
      index_index != std::string::npos) {
    return_string.replace(index_index, index_string.size(),
                          std::to_string(index));
  }

  if (auto value_index = return_string.find(value_string);
      value_index != std::string::npos) {
    std::ostringstream value_stream;
    value_stream << std::scientific << std::setprecision(precision_) << value;
    return_string.replace(value_index, value_string.size(), value_stream.str());
  }

  return return_string;
}

bool IntDoublePairToString::is_registered_ =
    ConverterIFactory<std::pair<int, double>, std::string>::get()
    .RegisterConstructor(
        converter::ConverterName::kIntDoublePairToString,
        [](){
          std::unique_ptr<ConverterI<std::pair<int, double>, std::string>>
              return_ptr = std::make_unique<IntDoublePairToString>();
          return return_ptr;
        });

} // namespace to_string

} // namespace converter

} // namespace instrumentation

} // namespace bart
