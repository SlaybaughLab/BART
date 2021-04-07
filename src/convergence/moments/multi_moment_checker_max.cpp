#include "convergence/moments/multi_moment_checker_max.hpp"

#include <stdexcept>

namespace bart::convergence::moments {

auto MultiMomentCheckerMax::IsConverged(const MomentsMap &current_iteration,
                                        const MomentsMap &previous_iteration) -> bool {
  AssertThrow(current_iteration.size() > 0, dealii::ExcMessage("Current iteration moments map is empty"));
  AssertThrow(previous_iteration.size() > 0, dealii::ExcMessage("Previous iteration moments map is empty"));
  AssertThrow(previous_iteration.size() == current_iteration.size(),
              dealii::ExcMessage("Current and previous iterations must be the same size"));
  is_converged_ = true;
  delta_ = std::nullopt;
  failed_index_ = std::nullopt;

  for (const auto &[index, previous_moment] : previous_iteration) {
    const auto &[group, harmonic_l, harmonic_m] = index;

    // Check that l = m = 0 (scalar flux)
    if (harmonic_l == 0 && harmonic_m == 0) {
      try {
        auto current_moment = current_iteration.at(index);

        if (!checker_->IsConverged(current_moment, previous_moment)) {
          is_converged_ = false;

          const double delta = checker_->delta().value_or(0);

          if (delta > delta_.value_or(0)) {
            delta_ = delta;
            failed_index_ = group;
          }
        }
      } catch (std::out_of_range &exc) {
        AssertThrow(false, dealii::ExcMessage("Current iteration lacks a group that previous iteration had"));
      }
    }
  }

  return is_converged_;
}

} // namespace bart::convergence::moments