#include "instrumentation/converter/int_double_pair_to_string.h"

namespace bart {

namespace instrumentation {

namespace converter {

namespace  {

std::string default_output_format{"${INDEX}, ${VALUE}\n"};
std::map<IntDoublePairToStringOutputTerm, std::string>
    default_output_term_to_string_map{
    {kValue, "${VALUE}"},
    {kIndex, "${INDEX}"}};
}

IntDoublePairToString::IntDoublePairToString()
    : ToStringConverter<std::pair<int, double>, OutputTerm>(
    default_output_format, default_output_term_to_string_map){}

std::string IntDoublePairToString::Convert(const std::pair<int,
                                                           double> &input) const {
  return std::string();
}

} // namespace converter

} // namespace instrumentation

} // namespace bart
