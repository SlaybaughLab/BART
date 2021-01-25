#ifndef BART_SRC_ITERATION_INITIALIZER_INITIALIZE_FIXED_TERMS_RESET_MOMENTS_HPP_
#define BART_SRC_ITERATION_INITIALIZER_INITIALIZE_FIXED_TERMS_RESET_MOMENTS_HPP_

#include "iteration/initializer/initialize_fixed_terms.h"

namespace bart::iteration::initializer {

class InitializeFixedTermsResetMoments : public InitializeFixedTerms {
 public:
  using FixedUpdater = InitializeFixedTerms::FixedUpdaterType;
  InitializeFixedTermsResetMoments(std::shared_ptr<FixedUpdater> fixed_updater_ptr, const int total_groups,
                                   const int total_angles)
      : InitializeFixedTerms(fixed_updater_ptr, total_groups, total_angles) {}

  void Initialize(system::System &system) override {
    for (auto& [index, moment] : *system.current_moments)
      moment = 1.0;
    for (auto& [index, moment] : *system.previous_moments)
      moment = 1.0;
    system.k_effective = std::nullopt;
    InitializeFixedTerms::Initialize(system);
  }

};

} // namespace bart::iteration::initializer

#endif //BART_SRC_ITERATION_INITIALIZER_INITIALIZE_FIXED_TERMS_RESET_MOMENTS_HPP_
