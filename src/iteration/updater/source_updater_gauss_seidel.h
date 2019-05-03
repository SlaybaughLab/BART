#ifndef BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_GAUSS_SEIDEL_H_
#define BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_GAUSS_SEIDEL_H_

#include "iteration/updater/source_updater.h"

namespace bart {

namespace iteration {

namespace updater {

class SourceUpdaterGaussSeidel : public SourceUpdater {
 public:
  void UpdateScatteringSource(System& system);
};

} // namespace updater

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_GAUSS_SEIDEL_H_