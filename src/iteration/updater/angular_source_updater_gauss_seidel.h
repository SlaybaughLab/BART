#ifndef BART_SRC_ITERATION_UPDATER_ANGULAR_SOURCE_UPDATER_GAUSS_SEIDEL_H_
#define BART_SRC_ITERATION_UPDATER_ANGULAR_SOURCE_UPDATER_GAUSS_SEIDEL_H_

#include <memory>

#include "system/system_types.h"
#include "iteration/updater/source_updater.h"

namespace bart {

namespace system {
struct System;
} // namespace system

namespace iteration {

namespace updater {

template <typename StamperType>
class AngularSourceUpdaterGaussSeidel : public SourceUpdater<StamperType> {
 public:

  void UpdateScatteringSource(system::System& system,
                              system::GroupNumber group,
                              system::AngleIndex angle) override {};
  void UpdateFissionSource(system::System& system,
                           system::GroupNumber group,
                           system::AngleIndex angle) override {};
};

} // namespace updater

} // namespace iteration

} // namespace bart

#endif //BART_SRC_ITERATION_UPDATER_ANGULAR_SOURCE_UPDATER_GAUSS_SEIDEL_H_
