#ifndef BART_SRC_FORMULATION_UPDATER_DIFFUSION_UPDATER_H_
#define BART_SRC_FORMULATION_UPDATER_DIFFUSION_UPDATER_H_

#include <memory>

#include "formulation/scalar/diffusion_i.h"
#include "formulation/stamper_i.h"
#include "formulation/updater/fixed_updater_i.h"
#include "formulation/updater/scattering_source_updater_i.h"
#include "formulation/updater/fission_source_updater_i.h"
#include "quadrature/quadrature_set_i.h"

namespace bart {

namespace formulation {

namespace updater {

template <int dim>
class DiffusionUpdater {
 public:
  using DiffusionFormulationType = formulation::scalar::DiffusionI<dim>;
  using StamperType = formulation::StamperI<dim>;
  DiffusionUpdater(std::unique_ptr<DiffusionFormulationType>,
                   std::unique_ptr<StamperType>);
  virtual ~DiffusionUpdater() = default;


  DiffusionFormulationType* formulation_ptr() const {
    return formulation_ptr_.get(); }
  StamperType* stamper_ptr() const { return stamper_ptr_.get(); }
 private:
  std::unique_ptr<DiffusionFormulationType> formulation_ptr_;
  std::unique_ptr<StamperType> stamper_ptr_;
};

} // namespace updater

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_UPDATER_DIFFUSION_UPDATER_H_
