#ifndef BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_GAUSS_SEIDEL_H_
#define BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_GAUSS_SEIDEL_H_

#include <memory>

#include "data/system.h"
#include "data/system/system_types.h"
#include "iteration/updater/source_updater.h"

namespace bart {

namespace iteration {

namespace updater {

/*! \brief Updates the source terms in a system using the provided stamper and
 *         the current iteration moments.
 *
 * Implementation of this class must be created for each stamper, as they may
 * have different methods for updating source terms. The entire system is passed
 * to the member functions, allowing for updating or accessing of any part of
 * the system.
 *
 * When moments are used to update the sources, the moments from the current
 * solve are used, making this a *Gauss-Seidel-like* update.
 *
 * @tparam StamperType
 */
template <typename StamperType>
class SourceUpdaterGaussSeidel : public SourceUpdater<StamperType> {
 public:

  using VariableTerms = data::system::VariableLinearTerms;

  explicit SourceUpdaterGaussSeidel(std::unique_ptr<StamperType> stamper_ptr)
      : SourceUpdater<StamperType>(std::move(stamper_ptr)) {};

  void UpdateScatteringSource(data::System& system,
                              data::system::GroupNumber group,
                              data::system::AngleIndex angle) override;
  void UpdateFissionSource(data::System& system,
                           data::system::GroupNumber group,
                           data::system::AngleIndex angle) override;
};

} // namespace updater

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_GAUSS_SEIDEL_H_