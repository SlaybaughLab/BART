#ifndef BART_SRC_ITERATION_INITIALIZER_SET_FIXED_TERMS_ONCE_H_
#define BART_SRC_ITERATION_INITIALIZER_SET_FIXED_TERMS_ONCE_H_

#include "iteration/initializer/set_fixed_terms.h"

namespace bart {

namespace iteration {

namespace initializer {

/*! \brief Initializes a system by setting the fixed terms once.
 *
 * Identical to iteration::initializer::SetFixedTerms but after the first
 * call to Initialize, subsequent calls will not do anything.
 *
 */
class SetFixedTermsOnce : public SetFixedTerms {
 public:
  SetFixedTermsOnce(std::unique_ptr<updater::FixedUpdaterI> fixed_updater_ptr,
                    const int total_groups,
                    const int total_angles)
      : SetFixedTerms(std::move(fixed_updater_ptr), total_groups, total_angles)
      {};


  void Initialize(system::System& sys) override {
    if (!initialize_was_called_) {
      SetFixedTerms::Initialize(sys);
      initialize_was_called_ = true;
    }
  }

  void set_initialize_was_called(const bool to_set) {
    initialize_was_called_ = to_set; }
  bool initialize_was_called() const { return initialize_was_called_; }
 protected:
  bool initialize_was_called_ = false;
};

} // namespace initializer

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_INITIALIZER_SET_FIXED_TERMS_ONCE_H_