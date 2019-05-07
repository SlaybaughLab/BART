#ifndef BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_GAUSS_SEIDEL_H_
#define BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_GAUSS_SEIDEL_H_

#include <memory>

#include "data/system.h"
#include "data/system/system_types.h"
#include "iteration/updater/source_updater.h"


namespace bart {

namespace iteration {

namespace updater {

template <typename StamperType>
class SourceUpdaterGaussSeidel : public SourceUpdater<StamperType> {
 public:
  explicit SourceUpdaterGaussSeidel(std::unique_ptr<StamperType> stamper_ptr)
      : SourceUpdater<StamperType>(std::move(stamper_ptr)) {};

  void UpdateScatteringSource(data::System& system,
                              data::system::Index index) override;
};

} // namespace updater

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_GAUSS_SEIDEL_H_