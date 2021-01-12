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
  // Extract the angular flux for this group and all quadrature angles
  std::map<quadrature::QuadraturePointIndex, std::shared_ptr<dealii::Vector<double>>> group_angular_flux;
  for (auto& [angular_flux_index, vector_ptr] : angular_flux_storage_map_) {
    auto& [solution_energy_group, solution_angle_index] = angular_flux_index;
    if (solution_energy_group == energy_group)
      group_angular_flux.insert({quadrature::QuadraturePointIndex(solution_angle_index.get()), vector_ptr});
  }
  // Get scalar flux
  auto scalar_flux = this->high_order_moments_->GetMoment({energy_group.get(), 0, 0});
  // Get current at all degrees of freedom
  auto current_at_global_dofs = this->integrated_flux_calculator_ptr()->NetCurrent(group_angular_flux);
  // Break the net current at each degree of freedom into dim arrays that represent each component at each dof
  std::array<Vector, dim> current_directional_components_at_global_dofs;
  for (int dir = 0; dir < dim; ++dir) {
    const int global_dofs = current_at_global_dofs.size();
    Vector current_component_at_dofs(global_dofs);
    for (int i = 0; i < global_dofs; ++i) {
      current_component_at_dofs[i] = current_at_global_dofs.at(i)[dir];
    }
    current_directional_components_at_global_dofs[dir] = current_component_at_dofs;
  }


  const auto drift_diffusion_term_function = [=, this](formulation::FullMatrix& cell_matrix,
                                                       const CellPtr& cell_ptr) -> void {
    drift_diffusion_formulation_ptr_->FillCellDriftDiffusionTerm(
        cell_matrix, cell_ptr, energy_group, scalar_flux, current_directional_components_at_global_dofs);
  };
  this->fixed_matrix_functions_.push_back(drift_diffusion_term_function);

  auto boundary_factory_function = [&](const dealii::Tensor<1, dim>& normal_tensor) -> dealii::Vector<double> {
    dealii::Vector<double> normal_vector(dim);
    for (int i = 0; i < dim; ++i) {
      normal_vector[i] = normal_tensor[i];
    }
    auto directional_current = integrated_flux_calculator_ptr_->DirectionalCurrent(group_angular_flux, normal_vector);
    auto directional_flux = integrated_flux_calculator_ptr_->DirectionalFlux(group_angular_flux, normal_vector);
    dealii::Vector<double> return_vector(directional_current.size());
    for (unsigned int i = 0; i < directional_current.size(); ++i) {
      if (directional_flux.at(i) == 0) {
        return_vector[i] = 0;
      } else {
        return_vector[i] = (directional_current.at(i)/directional_flux.at(i));
      }
    }
    return return_vector;
  };

  const auto drift_diffusion_boundary_function = [&, this](formulation::FullMatrix& cell_matrix,
                                                           const domain::FaceIndex face_index,
                                                           const CellPtr& cell_ptr) -> void {
    using BoundaryType = formulation::BoundaryType;
    problem::Boundary boundary = static_cast<problem::Boundary>(cell_ptr->face(face_index.get())->boundary_id());
    BoundaryType boundary_type{ BoundaryType::kVacuum };

    if (this->reflective_boundaries_.count(boundary) == 1)
      boundary_type = BoundaryType::kReflective;
    drift_diffusion_formulation_ptr_->FillCellBoundaryTerm(cell_matrix, cell_ptr, face_index, boundary_type,
                                                           boundary_factory_function);
  };
  this->fixed_matrix_boundary_functions_.push_back(drift_diffusion_boundary_function);
}

template class DriftDiffusionUpdater<1>;
template class DriftDiffusionUpdater<2>;
template class DriftDiffusionUpdater<3>;

} // namespace bart::formulation::updater
