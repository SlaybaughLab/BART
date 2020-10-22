#include "instrumentation/converter/convert_to_string/int_vector_complex_pair_to_string.h"

#include "instrumentation/converter/factory.hpp"

#include <iomanip>

namespace bart {

namespace instrumentation {

namespace converter {

namespace convert_to_string {

namespace  {
using OutputTerm = IntVectorComplexPairToStringOutputTerm;
std::string default_output_format{"${VECTOR}\n"};
std::map<OutputTerm, std::string> default_output_term_to_string_map{
    {OutputTerm::kInt,     "${INT}"},
    {OutputTerm::kVector, "${VECTOR}"},
};
using VectorTerm = IntVectorComplexPairToStringVectorOutputTerm;
std::string default_vector_entry_output_format{"(${REAL}+${IMAG}j)"};
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
  using VectorTerm = IntVectorComplexPairToStringVectorOutputTerm;
  std::string return_string = output_format_;
  const auto& [integer, vector] = input;

  std::string integer_string{output_term_to_string_map_.at(OutputTerm::kInt)};
  if (auto integer_index = return_string.find(integer_string);
      integer_index != std::string::npos) {
    return_string.replace(integer_index, integer_string.size(),
                          std::to_string(integer));
  }

  std::string vector_string{output_term_to_string_map_.at(OutputTerm::kVector)};
  if (auto vector_index = return_string.find(vector_string);
      vector_index != std::string::npos) {
    std::string full_vector_string{""};
    for (int i = 0; i < static_cast<int>(vector.size()); ++i) {
      std::string vector_entry_string{vector_output_format_};
      for (const auto& [term, string] : vector_term_to_string_map_) {
        if (auto term_index = vector_entry_string.find(string);
            term_index != std::string::npos) {
          std::ostringstream value_stream;
          if (term == VectorTerm::kReal || term == VectorTerm::kImag) {
            value_stream << std::scientific << std::setprecision(precision_);
          }
          switch (term) {
            case (VectorTerm::kReal): {
              value_stream << vector.at(i).real();
              break;
            }
            case (VectorTerm::kImag): {
              value_stream << vector.at(i).imag();
              break;
            }
            case (VectorTerm::kIndex):
            default: {
                  value_stream << i;
                  break;
            }
          }
          vector_entry_string.replace(term_index, string.size(), value_stream.str());
        }
      }
      full_vector_string += vector_entry_string;
      if ((i + 1) != static_cast<int>(vector.size()))
        full_vector_string += vector_delimiter_;
    }
    return_string.replace(vector_index, vector_string.size(), full_vector_string);
  }

  return return_string;
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

bool IntVectorComplexPairToString::is_registered_ = ConverterIFactory<InputType, std::string>::get()
    .RegisterConstructor(
        converter::ConverterName::kIntVectorComplexPairToString,
        []() {
          std::unique_ptr<ConverterI<InputType, std::string>> return_ptr =
              std::make_unique<IntVectorComplexPairToString>();
          return return_ptr;
        });


} // namespace convert_to_string

} // namespace converter

} // namespace instrumentation

} // namespace bart
