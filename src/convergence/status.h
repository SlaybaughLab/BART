#ifndef BART_SRC_CONVERGENCE_STATUS_H_
#define BART_SRC_CONVERGENCE_STATUS_H_

#include <optional>

namespace bart {

namespace convergence {

/*! Contains the status of a convergence check */
struct Status {
  int iteration_number = 0;
  int max_iterations = 0;
  bool is_complete = false;
  std::optional<int> failed_index = std::nullopt;
  std::optional<double> delta = std::nullopt;
};

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_STATUS_H_    
