#include "instrumentation/converter/int_vector_complex_pair_to_string.h"

namespace bart {

namespace instrumentation {

namespace converter {

namespace  {
using OutputTerm = IntVectorComplexPairToStringOutputTerm;
std::string default_output_format{"${INT}, (${REAL}, ${IMAG})"};
std::map<OutputTerm, std::string> default_output_term_to_string_map{
    {OutputTerm::kInt,         "${INT}"},
    {OutputTerm::kVectorIndex, "${INDEX}"},
    {OutputTerm::kReal,        "${REAL}"},
    {OutputTerm::kImag,        "${IMAG}"}
};

} // namespace


IntVectorComplexPairToString::IntVectorComplexPairToString()
    :  ToStringConverter<InputType, OutputTerm>(
        default_output_format, default_output_term_to_string_map) {}

std::string IntVectorComplexPairToString::Convert(
    const std::pair<int,std::vector<std::complex<double>>> &input) const {
  return std::string();
}

} // namespace converter

} // namespace instrumentation

} // namespace bart
