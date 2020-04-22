#ifndef BART_SRC_FORMULATION_UPDATER_TESTS_SAAF_UPDATER_H_
#define BART_SRC_FORMULATION_UPDATER_TESTS_SAAF_UPDATER_H_

#include <memory>

#include "formulation/angular/self_adjoint_angular_flux_i.h"
#include "formulation/stamper_i.h"
#include "formulation/updater/fixed_updater_i.h"
#include "formulation/updater/scattering_source_updater_i.h"
#include "formulation/updater/fission_source_updater_i.h"
#include "quadrature/quadrature_set_i.h"

namespace bart {

namespace formulation {

namespace updater {

template <int dim>
class SAAFUpdater :
    public FixedUpdaterI,
    public ScatteringSourceUpdaterI,
    public FissionSourceUpdaterI {
 public:
  using SAAFFormulationType = formulation::angular::SelfAdjointAngularFluxI<dim>;
  using StamperType = formulation::StamperI<dim>;
  using QuadratureSetType = quadrature::QuadratureSetI<dim>;
  SAAFUpdater(std::unique_ptr<SAAFFormulationType>,
              std::unique_ptr<StamperType>,
              const std::shared_ptr<QuadratureSetType>&);

  void UpdateFixedTerms(system::System &to_update,
                        system::EnergyGroup group,
                        quadrature::QuadraturePointIndex index) override;
  void UpdateFissionSource(system::System &to_update,
                           system::EnergyGroup group,
                           quadrature::QuadraturePointIndex index) override;
  void UpdateScatteringSource(system::System &to_update,
                              system::EnergyGroup group,
                              quadrature::QuadraturePointIndex index) override;

  SAAFFormulationType* formulation_ptr() const {return formulation_ptr_.get();};
  StamperType* stamper_ptr() const {return stamper_ptr_.get();};
  QuadratureSetType* quadrature_set_ptr() const {
    return quadrature_set_ptr_.get();};
 private:
  std::unique_ptr<SAAFFormulationType> formulation_ptr_;
  std::unique_ptr<StamperType> stamper_ptr_;
  std::shared_ptr<QuadratureSetType> quadrature_set_ptr_;
};

} // namespace updater

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_UPDATER_TESTS_SAAF_UPDATER_H_
