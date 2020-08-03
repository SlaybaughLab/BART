#include "instrumentation/converter/convergence_to_string.h"

#include <sstream>

namespace bart {

namespace instrumentation {

namespace converter {

ConvergenceToString::ConvergenceToString()
    : ToStringConverter<convergence::Status, ConvergenceToStringOutputTerm>(
    "Iteration: ${ITERATION_NUM}/${ITERATION_MAX}, delta: ${DELTA}, index: ${INDEX}\n",
    {{kIterationNum, "${ITERATION_NUM}"},
     {kIterationMax, "${ITERATION_MAX}"},
     {kDelta, "${DELTA}"},
     {kIndex, "${INDEX}"}}) {}

std::string ConvergenceToString::Convert(const convergence::Status &to_convert) const {
  auto return_string = output_format_;

  std::string delta_string{null_character_}, index_string{null_character_};
  if (to_convert.delta.has_value()) {
    std::ostringstream delta_stream;
    delta_stream << to_convert.delta.value();
    delta_string = delta_stream.str();
  }
  if (to_convert.failed_index.has_value()) {
    index_string = std::to_string(to_convert.failed_index.value());
  }


  std::map<OutputTerm, std::string> output_term_string_map{
      {OutputTerm::kIterationNum, std::to_string(to_convert.iteration_number)},
      {OutputTerm::kIterationMax, std::to_string(to_convert.max_iterations)},
      {OutputTerm::kIndex, index_string},
      {OutputTerm::kDelta, delta_string}
  };

  for (const auto& [term, value] : output_term_string_map) {
    std::string string_to_find = output_term_to_string_map_.at(term);
    if (auto index = return_string.find(string_to_find);
        index != std::string::npos) {
      return_string.replace(index, string_to_find.size(), value);
    }
  }

  return return_string;
}

std::string ConvergenceToString::SetOutputFormat(
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
}

} // namespace converter

} // namespace instrumentation

} // namespace bart
