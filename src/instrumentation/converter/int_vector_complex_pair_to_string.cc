#include "instrumentation/converter/int_vector_complex_pair_to_string.h"

namespace bart {

namespace instrumentation {

namespace converter {

namespace  {
using OutputTerm = IntVectorComplexPairToStringOutputTerm;
std::string default_output_format{"${INT}, ${VECTOR}\n"};
std::map<OutputTerm, std::string> default_output_term_to_string_map{
    {OutputTerm::kInt,     "${INT}"},
    {OutputTerm::kVector, "${VECTOR}"},
};
using VectorTerm = IntVectorComplexPairToStringVectorOutputTerm;
std::string default_vector_entry_output_format{"(${REAL}, ${IMAG}$)"};
std::map<VectorTerm, std::string> default_vector_entry_output_to_string_map{
    {VectorTerm::kIndex, "${INDEX}"},
    {VectorTerm::kReal, "${REAL}"},
    {VectorTerm::kImag, "${IMAG}"}
};

} // namespace

IntVectorComplexPairToString::IntVectorComplexPairToString()
    : vector_output_format_(default_vector_entry_output_format),
      vector_term_to_string_map_(default_vector_entry_output_to_string_map),
      ToStringConverter<InputType, OutputTerm>(
          default_output_format, default_output_term_to_string_map) {}

std::string IntVectorComplexPairToString::Convert(
    const std::pair<int,std::vector<std::complex<double>>> &input) const {
  return std::string();
}
std::string IntVectorComplexPairToString::SetVectorOutputFormat(
    std::vector<std::variant<VectorTerm,std::string>> output_format_parts) {

  std::ostringstream new_format_stream;
  for (const auto output_part : output_format_parts) {
    try {
      new_format_stream << vector_term_to_string_map_.at(
          std::get<VectorTerm>(output_part));
    } catch (const std::bad_variant_access&) {
      new_format_stream << std::get<std::string>(output_part);
    }
  }
  vector_output_format_ = new_format_stream.str();
  return vector_output_format_;
}

} // namespace converter

} // namespace instrumentation

} // namespace bart
