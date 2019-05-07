#ifndef BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_I_H_
#define BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_I_H_

#include "data/system/system_types.h"
#include "data/system/system.h"

namespace bart {

namespace iteration {

namespace updater {

class SourceUpdaterI {
 public:
  using System = data::system::System;
  virtual ~SourceUpdaterI() = default;

  virtual void UpdateScatteringSource(System& system,
                                      data::system::Index index) = 0;
};

} // namespace updater

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_I_H_