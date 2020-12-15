#include "drift_diffusion_updater.hpp"

namespace bart::formulation::updater {

template<int dim>
DriftDiffusionUpdater<dim>::DriftDiffusionUpdater(
    std::unique_ptr<DiffusionFormulation> diffusion_formulation_ptr,
    std::unique_ptr<DriftDiffusionFormulation> drift_diffusion_formulation_ptr,
    std::shared_ptr<Stamper> stamper_ptr,
    std::unique_ptr<IntegratedFluxCalculator> integrated_flux_calculator_ptr,
    AngularFluxStorageMap& angular_flux_storage_map,
    std::unordered_set<problem::Boundary> reflective_boundaries)
    : DiffusionUpdater<dim>(std::move(diffusion_formulation_ptr), stamper_ptr, reflective_boundaries),
        angular_flux_storage_map_(angular_flux_storage_map),
        drift_diffusion_formulation_ptr_(std::move(drift_diffusion_formulation_ptr)),
        integrated_flux_calculator_ptr_(std::move(integrated_flux_calculator_ptr)) {
  std::string call_location{ "DriftDiffusionUpdaterConstructor"};
  AssertPointerNotNull(drift_diffusion_formulation_ptr_.get(), "drift diffusion formulation", call_location);
  AssertPointerNotNull(integrated_flux_calculator_ptr_.get(), "integrated flux calculator", call_location);
}
template<int dim>
auto DriftDiffusionUpdater<dim>::SetUpFixedFunctions(system::System& system,
                                                     system::EnergyGroup energy_group,
                                                     quadrature::QuadraturePointIndex quadrature_point_index) -> void {
  DiffusionUpdater<dim>::SetUpFixedFunctions(system, energy_group, quadrature_point_index);
}

template class DriftDiffusionUpdater<1>;
template class DriftDiffusionUpdater<2>;
template class DriftDiffusionUpdater<3>;

} // namespace bart::formulation::updater
