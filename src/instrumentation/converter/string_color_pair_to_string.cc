#include "instrumentation/converter/string_color_pair_to_string.h"

namespace bart {

namespace instrumentation {

namespace converter {

namespace  {
using OutputTerm = StringColorPairToStringOutputTerm;
std::string default_output_format{"${COLOR_CODE}${STRING}${COLOR_RESET_CODE}"};
std::map<StringColorPairToStringOutputTerm, std::string>
    default_output_term_to_string_map{
    {OutputTerm::kColorCode,  "${COLOR_CODE}"},
    {OutputTerm::kString,     "${STRING}"},
    {OutputTerm::kColorReset, "${COLOR_RESET_CODE}"}};

} // namespace

StringColorPairToString::StringColorPairToString()
: ToStringConverter<std::pair<std::string, utility::Color>, OutputTerm>(
    default_output_format, default_output_term_to_string_map) {}

std::string StringColorPairToString::Convert(const std::pair<std::string,
                                                             utility::Color> &input) const {
  return std::string();
}

} // namespace converter

} // namespace instrumentation

} // namespace bart
