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
  SetFixedTerms(std::unique_ptr<updater::FixedUpdaterI> fixed_updater_ptr)
      : fixed_updater_ptr_(std::move(fixed_updater_ptr)) {};

  virtual ~SetFixedTerms() = default;
  virtual void Initialize(data::System& sys) override {};

  updater::FixedUpdaterI* GetUpdater() const {
    return fixed_updater_ptr_.get(); };

 protected:
  std::unique_ptr<updater::FixedUpdaterI> fixed_updater_ptr_;
};

} // namespace initializer

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_INITIALIZER_SET_FIXED_TERMS_H_