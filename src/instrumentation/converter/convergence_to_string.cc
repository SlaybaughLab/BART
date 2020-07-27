#include "instrumentation/converter/convergence_to_string.h"

#include <sstream>

namespace bart {

namespace instrumentation {

namespace converter {

std::string ConvergenceToString::Convert(const convergence::Status &to_convert) const {
  auto return_string = output_format_;
  // Iteration number
  std::string iteration_num_string = output_term_to_string_map_.at(OutputTerm::kIterationNum);
  return_string.replace(return_string.find(iteration_num_string),
                 iteration_num_string.size(),
                 std::to_string(to_convert.iteration_number));
  // Max iterations
  std::string iteration_max_string = output_term_to_string_map_.at(OutputTerm::kIterationMax);
  return_string.replace(return_string.find(iteration_max_string),
                 iteration_max_string.size(),
                 std::to_string(to_convert.max_iterations));
  // Delta
  std::ostringstream delta_stream;
  delta_stream << to_convert.delta.value();
  std::string delta_string = output_term_to_string_map_.at(OutputTerm::kDelta);
  return_string.replace(return_string.find(delta_string),
                 delta_string.size(),
                 delta_stream.str());
  // Index
  std::string index_string = output_term_to_string_map_.at(OutputTerm::kIndex);
  return_string.replace(return_string.find(index_string),
                 index_string.size(),
                 std::to_string(to_convert.failed_index.value()));

  return return_string;
}
} // namespace converter

} // namespace instrumentation

} // namespace bart
