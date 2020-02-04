#ifndef BART_SRC_ITERATION_UPDATER_ANGULAR_SOURCE_UPDATER_GAUSS_SEIDEL_H_
#define BART_SRC_ITERATION_UPDATER_ANGULAR_SOURCE_UPDATER_GAUSS_SEIDEL_H_

#include <memory>

#include "system/system_types.h"
#include "iteration/updater/source_updater.h"
#include "quadrature/quadrature_set_i.h"

namespace bart {

namespace system {
struct System;
} // namespace system

namespace iteration {

namespace updater {

template <typename StamperType>
class AngularSourceUpdaterGaussSeidel : public SourceUpdater<StamperType> {
 public:

  static constexpr int dim = StamperType::dimension;
  using QuadratureSetType = quadrature::QuadratureSetI<dim>;

  AngularSourceUpdaterGaussSeidel(
      std::shared_ptr<StamperType> stamper_ptr,
      std::shared_ptr<QuadratureSetType>);

  void UpdateScatteringSource(system::System& system,
                              system::GroupNumber group,
                              system::AngleIndex angle) override {};
  void UpdateFissionSource(system::System& system,
                           system::GroupNumber group,
                           system::AngleIndex angle) override {};

  StamperType* stamper_ptr() const { return  this->stamper_ptr_.get(); }
  QuadratureSetType* quadrature_set_ptr() const {
    return quadrature_set_ptr_.get(); }

 private:
  std::shared_ptr<QuadratureSetType> quadrature_set_ptr_ = nullptr;
};

} // namespace updater

} // namespace iteration

} // namespace bart

#endif //BART_SRC_ITERATION_UPDATER_ANGULAR_SOURCE_UPDATER_GAUSS_SEIDEL_H_
