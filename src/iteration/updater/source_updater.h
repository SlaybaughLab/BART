#ifndef BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_H_
#define BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_H_

#include "iteration/updater/source_updater_i.h"

namespace bart {

namespace iteration {

namespace updater {

class SourceUpdater : public SourceUpdaterI {
 public:
  virtual ~SourceUpdater() = default;

};

} // namespace updater

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_H_