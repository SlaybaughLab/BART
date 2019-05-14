#ifndef BART_SRC_ITERATION_UPDATER_FIXED_UPDATER_H_
#define BART_SRC_ITERATION_UPDATER_FIXED_UPDATER_H_

#include "iteration/updater/fixed_updater_i.h"

namespace bart {

namespace iteration {

namespace updater {

class FixedUpdater : public FixedUpdaterI {
 public:
  virtual ~FixedUpdater() = default;

  void UpdateFixedTerms(data::System& system,
                        data::system::GroupNumber group,
                        data::system::AngleIndex angle) override {};
};

} // namespace updater

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_UPDATER_FIXED_UPDATER_H_