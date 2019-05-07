#ifndef BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_I_H_
#define BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_I_H_

#include "data/system/system_types.h"
#include "data/system.h"

namespace bart {

namespace iteration {

namespace updater {

class SourceUpdaterI {
 public:
  virtual ~SourceUpdaterI() = default;

  virtual void UpdateScatteringSource(data::System& system,
                                      data::system::GroupNumber group,
                                      data::system::AngleIndex angle) = 0;
};

} // namespace updater

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_I_H_