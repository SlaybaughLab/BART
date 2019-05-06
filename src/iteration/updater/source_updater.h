#ifndef BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_H_
#define BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_H_

#include <memory>

#include "iteration/updater/source_updater_i.h"
#include "formulation/cfem_stamper_i.h"

namespace bart {

namespace iteration {

namespace updater {

template <typename StamperType>
class SourceUpdater : public SourceUpdaterI {
 public:
  explicit SourceUpdater(std::unique_ptr<StamperType> stamper_ptr)
      : stamper_ptr_(std::move(stamper_ptr)) {};
  virtual ~SourceUpdater() override = default;

 protected:
  std::unique_ptr<StamperType> stamper_ptr_;

};

} // namespace updater

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_H_