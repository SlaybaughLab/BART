#ifndef BART_SRC_FORMULATION_UPDATER_DRIFT_DIFFUSION_UPDATER_HPP_
#define BART_SRC_FORMULATION_UPDATER_DRIFT_DIFFUSION_UPDATER_HPP_

#include "diffusion_updater.hpp"
#include "formulation/scalar/drift_diffusion_i.hpp"
#include "formulation/stamper_i.h"
#include "quadrature/calculators/drift_diffusion_integrated_flux_i.hpp"
#include "system/solution/solution_types.h"
#include "utility/has_dependencies.h"

namespace bart::formulation::updater {

template <int dim>
class DriftDiffusionUpdater : public DiffusionUpdater<dim>, public utility::HasDependencies {
 public:
  using AngularFluxStorageMap = system::solution::EnergyGroupToAngularSolutionPtrMap;
  using DiffusionFormulation = typename DiffusionUpdater<dim>::DiffusionFormulationType;
  using DriftDiffusionFormulation = formulation::scalar::DriftDiffusionI<dim>;
  using IntegratedFluxCalculator = quadrature::calculators::DriftDiffusionIntegratedFluxI;
  using Stamper = typename DiffusionUpdater<dim>::StamperType;

  DriftDiffusionUpdater(std::unique_ptr<DiffusionFormulation>,
                        std::unique_ptr<DriftDiffusionFormulation>,
                        std::shared_ptr<Stamper>,
                        std::unique_ptr<IntegratedFluxCalculator>,
                        AngularFluxStorageMap&,
                        std::unordered_set<problem::Boundary> reflective_boundaries = {});
  virtual ~DriftDiffusionUpdater() = default;

  auto angular_flux_storage_map() const -> AngularFluxStorageMap {
    return angular_flux_storage_map_; }
  auto drift_diffusion_formulation_ptr() const -> DriftDiffusionFormulation* {
    return drift_diffusion_formulation_ptr_.get(); }
  auto integrated_flux_calculator_ptr() const -> IntegratedFluxCalculator* {
    return integrated_flux_calculator_ptr_.get(); }
 protected:
  AngularFluxStorageMap angular_flux_storage_map_{};
  std::unique_ptr<DriftDiffusionFormulation> drift_diffusion_formulation_ptr_{ nullptr };
  std::unique_ptr<IntegratedFluxCalculator> integrated_flux_calculator_ptr_{ nullptr };
};

} // namespace bart::formulation::updater

#endif //BART_SRC_FORMULATION_UPDATER_DRIFT_DIFFUSION_UPDATER_HPP_
