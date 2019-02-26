#include "multi_checker.h"

namespace bart {

namespace convergence {

namespace flux {

std::optional<int> MultiChecker::GetFailedIndex() const {
  if (is_converged_)
    return {};
  else
    return failed_index_;
}

} // namespace flux

} // namespace convergence

} // namespace bart
