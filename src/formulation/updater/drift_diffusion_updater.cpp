#include "drift_diffusion_updater.hpp"

namespace bart::formulation::updater {

template<int dim>
DriftDiffusionUpdater<dim>::DriftDiffusionUpdater(
    std::unique_ptr<DiffusionFormulation> diffusion_formulation_ptr,
    std::unique_ptr<DriftDiffusionFormulation> drift_diffusion_formulation_ptr,
    std::shared_ptr<Stamper> stamper_ptr,
    std::unique_ptr<IntegratedFluxCalculator> integrated_flux_calculator_ptr,
    std::shared_ptr<HighOrderMoments> high_order_moments,
    AngularFluxStorageMap& angular_flux_storage_map,
    std::unordered_set<problem::Boundary> reflective_boundaries)
    : DiffusionUpdater<dim>(std::move(diffusion_formulation_ptr), stamper_ptr, reflective_boundaries),
        angular_flux_storage_map_(angular_flux_storage_map),
        high_order_moments_(high_order_moments),
        drift_diffusion_formulation_ptr_(std::move(drift_diffusion_formulation_ptr)),
        integrated_flux_calculator_ptr_(std::move(integrated_flux_calculator_ptr)) {
  std::string call_location{ "DriftDiffusionUpdaterConstructor"};
  AssertPointerNotNull(drift_diffusion_formulation_ptr_.get(), "drift diffusion formulation", call_location);
  AssertPointerNotNull(integrated_flux_calculator_ptr_.get(), "integrated flux calculator", call_location);
  AssertPointerNotNull(high_order_moments_.get(), "higher order moments", call_location);
}
template<int dim>
auto DriftDiffusionUpdater<dim>::SetUpFixedFunctions(system::System& system,
                                                     system::EnergyGroup energy_group,
                                                     quadrature::QuadraturePointIndex quadrature_point_index) -> void {
  DiffusionUpdater<dim>::SetUpFixedFunctions(system, energy_group, quadrature_point_index);
  using CellPtr = domain::CellPtr<dim>;
  std::map<quadrature::QuadraturePointIndex, std::shared_ptr<dealii::Vector<double>>> group_angular_flux;
  for (auto& [angular_flux_index, vector_ptr] : angular_flux_storage_map_) {
    auto& [solution_energy_group, solution_angle_index] = angular_flux_index;
    if (solution_energy_group == energy_group)
      group_angular_flux.insert({quadrature::QuadraturePointIndex(solution_angle_index.get()), vector_ptr});
  }
  auto integrated_angular_flux = integrated_flux_calculator_ptr_->Integrate(group_angular_flux);
  auto scalar_flux = this->high_order_moments_->GetMoment({energy_group.get(), 0, 0});

  const auto drift_diffusion_term_function = [=, this](formulation::FullMatrix& cell_matrix,
                                                       const CellPtr& cell_ptr) -> void {
    drift_diffusion_formulation_ptr_->FillCellDriftDiffusionTerm(
        cell_matrix, cell_ptr, energy_group, scalar_flux, integrated_angular_flux);
  };
  this->fixed_matrix_functions_.push_back(drift_diffusion_term_function);
}

template class DriftDiffusionUpdater<1>;
template class DriftDiffusionUpdater<2>;
template class DriftDiffusionUpdater<3>;

} // namespace bart::formulation::updater
