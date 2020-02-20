#ifndef BART_SRC_FORMULATION_UPDATER_TESTS_SAAF_UPDATER_H_
#define BART_SRC_FORMULATION_UPDATER_TESTS_SAAF_UPDATER_H_

#include <memory>

#include "formulation/angular/cfem_self_adjoint_angular_flux_i.h"
#include "formulation/stamper_i.h"
#include "formulation/updater/fixed_updater_i.h"

namespace bart {

namespace formulation {

namespace updater {

template <int dim>
class SAAFUpdater {
 public:
  using SAAFFormulationType = formulation::angular::CFEMSelfAdjointAngularFluxI<dim>;
  using StamperType = formulation::StamperI<dim>;
  SAAFUpdater(std::unique_ptr<SAAFFormulationType>,
              std::unique_ptr<StamperType>);

  SAAFFormulationType* formulation_ptr() const {return formulation_ptr_.get();};
  StamperType* stamper_ptr() const {return stamper_ptr_.get();};
 private:
  std::unique_ptr<SAAFFormulationType> formulation_ptr_;
  std::unique_ptr<StamperType> stamper_ptr_;
};

} // namespace updater

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_UPDATER_TESTS_SAAF_UPDATER_H_
