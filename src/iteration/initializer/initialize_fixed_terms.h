#ifndef BART_SRC_ITERATION_INITIALIZER_INITIALIZE_FIXED_TERMS_H_
#define BART_SRC_ITERATION_INITIALIZER_INITIALIZE_FIXED_TERMS_H_

#include "formulation/updater/fixed_updater_i.h"
#include "iteration/initializer/initializer_i.h"

namespace bart {

namespace iteration {

namespace initializer {

class InitializeFixedTerms : public InitializerI {
 public:
  using FixedUpdaterType = formulation::updater::FixedUpdaterI;

  InitializeFixedTerms(std::unique_ptr<FixedUpdaterType> fixed_updater_ptr,
      int total_groups, int total_angles);
  void Initialize(system::System &sys) override;
  virtual ~InitializeFixedTerms() = default;

  FixedUpdaterType* fixed_updater_ptr() { return fixed_updater_ptr_.get(); };
  int total_groups() { return total_groups_; };
  int total_angles() { return total_angles_; };
 private:
  std::unique_ptr<FixedUpdaterType> fixed_updater_ptr_ = nullptr;
  const int total_groups_, total_angles_;
};

} // namespace initializer

} // namespace iteration

} // namespace bart

#endif //BART_SRC_ITERATION_INITIALIZER_INITIALIZE_FIXED_TERMS_H_
