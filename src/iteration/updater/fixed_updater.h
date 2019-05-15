#ifndef BART_SRC_ITERATION_UPDATER_FIXED_UPDATER_H_
#define BART_SRC_ITERATION_UPDATER_FIXED_UPDATER_H_

#include <memory>

#include "iteration/updater/fixed_updater_i.h"

namespace bart {

namespace iteration {

namespace updater {

template <typename StamperType>
class FixedUpdater : public FixedUpdaterI {
 public:
  explicit FixedUpdater(std::unique_ptr<StamperType> stamper_ptr);
  virtual ~FixedUpdater() = default;

  void UpdateFixedTerms(data::System& system,
                        data::system::GroupNumber group,
                        data::system::AngleIndex angle) override {};
  StamperType* GetStamperPtr() const {
    return stamper_ptr_.get();
  };

 protected:
  std::unique_ptr<StamperType> stamper_ptr_;
};

} // namespace updater

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_UPDATER_FIXED_UPDATER_H_