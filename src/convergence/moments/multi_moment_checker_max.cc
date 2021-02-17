#include "convergence/moments/multi_moment_checker_max.h"

#include <stdexcept>

namespace bart {

namespace convergence {

namespace moments {

bool MultiMomentCheckerMax::IsConverged(
    const system::moments::MomentsMap &current_iteration,
    const system::moments::MomentsMap &previous_iteration) {
  AssertThrow(current_iteration.size() > 0,
              dealii::ExcMessage("Current iteration moments map is empty"));
  AssertThrow(previous_iteration.size() > 0,
              dealii::ExcMessage("Previous iteration moments map is empty"));
  AssertThrow(previous_iteration.size() == current_iteration.size(),
              dealii::ExcMessage("Current and previous iterations must be the"
                                 "same size"));
  is_converged_ = true;

  for (auto &previous_pair : previous_iteration) {
    auto &[index, previous_moment] = previous_pair;

    // Check that l = m = 0 (scalar flux)
    if (index[1] == 0 && index[2] == 0) {
      try {
        auto current_moment = current_iteration.at(index);

        if (!checker_->IsConverged(current_moment, previous_moment)) {
          is_converged_ = false;

          double delta = checker_->delta().value_or(0);

          if (delta > delta_.value_or(0)) {
            delta_ = delta;
            failed_index_ = index[0];
          }
        }
      } catch (std::out_of_range &exc) {
        AssertThrow(false,
            dealii::ExcMessage("Current iteration lacks a group that previous"
                               "iteration had"));
      }
    }
  }

  if (is_converged_) {
    delta_ = std::nullopt;
    failed_index_ = std::nullopt;
  }

  return is_converged_;
}



} // namespace moments

} // namespace convergence

} // namespace bart