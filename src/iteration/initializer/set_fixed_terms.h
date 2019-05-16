#ifndef BART_SRC_ITERATION_INITIALIZER_SET_FIXED_TERMS_H_
#define BART_SRC_ITERATION_INITIALIZER_SET_FIXED_TERMS_H_

#include <memory>

#include "data/system.h"
#include "iteration/initializer/initializer_i.h"

namespace bart {

//namespace data { struct System; } // namespace data

namespace iteration {

namespace updater { class FixedUpdaterI;} // namespace updater

namespace initializer {
/*! \brief Initializes a system by setting the fixed terms.
 *
 * This class takes a class of type iteration::updater::FixedUpdaterI and
 * uses the UpdateFixedTerms method to initialize the system.
 *
 */
class SetFixedTerms : public InitializerI {
 public:
  SetFixedTerms(std::unique_ptr<updater::FixedUpdaterI> fixed_updater_ptr,
                const int total_groups,
                const int total_angles)
      : fixed_updater_ptr_(std::move(fixed_updater_ptr)),
        total_groups_(0),
        total_angles_(0) {};

  virtual ~SetFixedTerms() = default;
  virtual void Initialize(data::System& sys) override {};

  updater::FixedUpdaterI* fixed_updater_ptr() const {
    return fixed_updater_ptr_.get(); };

  int total_groups() const { return total_groups_; }
  int total_angles() const { return total_angles_; }

 protected:
  std::unique_ptr<updater::FixedUpdaterI> fixed_updater_ptr_;
  const int total_groups_;
  const int total_angles_;
};

} // namespace initializer

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_INITIALIZER_SET_FIXED_TERMS_H_