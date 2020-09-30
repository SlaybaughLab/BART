#include "instrumentation/converter/to_string/convergence_to_string.h"

#include <sstream>

#include "instrumentation/converter/factory.h"

namespace bart {

namespace instrumentation {

namespace converter {

namespace to_string {

namespace  {

using OutputTerm = ConvergenceToStringOutputTerm;

std::string default_output_format{"Iteration: ${ITERATION_NUM}/${ITERATION_MAX}, delta: ${DELTA}, index: ${INDEX}\n"};
std::map<ConvergenceToStringOutputTerm, std::string>
    default_output_term_to_string_map{{OutputTerm::kIterationNum, "${ITERATION_NUM}"},
                                      {OutputTerm::kIterationMax, "${ITERATION_MAX}"},
                                      {OutputTerm::kDelta, "${DELTA}"},
                                      {OutputTerm::kIndex, "${INDEX}"}};

} // namespace

ConvergenceToString::ConvergenceToString()
    : ToStringConverter<convergence::Status, ConvergenceToStringOutputTerm>(
    default_output_format, default_output_term_to_string_map) {}

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

bool ConvergenceToString::is_registered_ =
    ConverterIFactory<convergence::Status, std::string>::get()
    .RegisterConstructor(converter::ConverterName::kConvergenceToString,
        [](){
          std::unique_ptr<ConverterI<convergence::Status, std::string>>
              return_ptr = std::make_unique<ConvergenceToString>();
          return return_ptr;
        });
} // namespace to_string

} // namespace converter

} // namespace instrumentation

} // namespace bart
