#ifndef BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_GAUSS_SEIDEL_H_
#define BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_GAUSS_SEIDEL_H_

#include <memory>

#include "iteration/updater/source_updater.h"


namespace bart {

namespace iteration {

namespace updater {

template <typename StamperType>
class SourceUpdaterGaussSeidel : public SourceUpdater<StamperType> {
 public:
  using typename SourceUpdaterI::System;

  explicit SourceUpdaterGaussSeidel(std::unique_ptr<StamperType> stamper_ptr)
      : SourceUpdater<StamperType>(std::move(stamper_ptr)) {};
  void UpdateScatteringSource(System& system) override;
};

} // namespace updater

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_GAUSS_SEIDEL_H_