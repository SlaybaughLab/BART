#include "formulation/updater/saaf_updater.h"

namespace bart {

namespace formulation {

namespace updater {

template<int dim>
SAAFUpdater<dim>::SAAFUpdater(
    std::unique_ptr<SAAFFormulationType> formulation_ptr,
    std::unique_ptr<StamperType> stamper_ptr,
    const std::shared_ptr<QuadratureSetType>& quadrature_set_ptr)
    : formulation_ptr_(std::move(formulation_ptr)),
      stamper_ptr_(std::move(stamper_ptr)),
      quadrature_set_ptr_(quadrature_set_ptr) {
  AssertThrow(formulation_ptr_ != nullptr,
      dealii::ExcMessage("Error in constructor of SAAFUpdater, formulation "
                         "pointer passed is null"))
  AssertThrow(stamper_ptr_ != nullptr,
              dealii::ExcMessage("Error in constructor of SAAFUpdater, stamper "
                                 "pointer passed is null"))
  AssertThrow(quadrature_set_ptr_ != nullptr,
              dealii::ExcMessage("Error in constructor of SAAFUpdater, "
                                 "quadrature set pointer passed is null"))
  this->set_description("Self-adjoint angular flux updater",
                        utility::DefaultImplementation(true));
}

template<int dim>
SAAFUpdater<dim>::SAAFUpdater(
    std::unique_ptr<SAAFFormulationType> formulation_ptr,
    std::unique_ptr<StamperType> stamper_ptr,
    const std::shared_ptr<QuadratureSetType>& quadrature_set_ptr,
    const EnergyGroupToAngularSolutionPtrMap& angular_solution_ptr_map,
    const std::unordered_set<Boundary> reflective_boundaries)
    : SAAFUpdater(std::move(formulation_ptr), std::move(stamper_ptr),
                  quadrature_set_ptr) {
  reflective_boundaries_ = reflective_boundaries;
  angular_solution_ptr_map_ = angular_solution_ptr_map;
  this->set_description("Self-adjoint angular flux updater with reflective "
                        "boundaries",
                        utility::DefaultImplementation(true));
}

template<int dim>
void SAAFUpdater<dim>::UpdateBoundaryConditions(
    system::System &to_update,
    system::EnergyGroup group,
    quadrature::QuadraturePointIndex index) {
  using system::terms::VariableLinearTerms;
  auto boundary_vector_ptr = to_update.right_hand_side_ptr_->GetVariableTermPtr(
          {group.get(), index.get()},
          VariableLinearTerms::kReflectiveBoundaryCondition);
  const auto quadrature_point_ptr = quadrature_set_ptr_->GetQuadraturePoint(index);
  const auto reflected_quadrature_point_index =
      quadrature_set_ptr_->GetReflectionIndex(quadrature_point_ptr);

  AssertThrow(reflected_quadrature_point_index.has_value(),
      dealii::ExcMessage("Error in UpdateBoundaryConditions, passed "
                         "quadrature point has no reflection"))

  const auto incoming_flux = angular_solution_ptr_map_.at(
      system::SolutionIndex(group, reflected_quadrature_point_index.value()));
  if (incoming_flux->size() > 0) {
    auto reflective_boundary_term_function =
        [&](formulation::Vector &cell_vector,
            const domain::FaceIndex face_index,
            const domain::CellPtr <dim> &cell_ptr) -> void {
          if (IsOnReflectiveBoundary(cell_ptr, face_index)) {
            formulation_ptr_->FillReflectiveBoundaryLinearTerm(
                cell_vector,
                cell_ptr,
                face_index,
                quadrature_point_ptr,
                *incoming_flux);
          }
        };
    *boundary_vector_ptr = 0;
    stamper_ptr_->StampBoundaryVector(*boundary_vector_ptr,
                                      reflective_boundary_term_function);
  }
}

template<int dim>
void SAAFUpdater<dim>::UpdateFixedTerms(
    system::System &to_update,
    system::EnergyGroup group,
    quadrature::QuadraturePointIndex index) {
  auto fixed_matrix_ptr =
      to_update.left_hand_side_ptr_->GetFixedTermPtr({group.get(), index.get()});
  auto quadrature_point_ptr = quadrature_set_ptr_->GetQuadraturePoint(index);
  auto streaming_term_function =
      [&](formulation::FullMatrix& cell_matrix,
          const domain::CellPtr<dim>& cell_ptr) -> void {
        formulation_ptr_->FillCellStreamingTerm(cell_matrix, cell_ptr,
                                                quadrature_point_ptr, group);
      };
  auto collision_term_function =
      [&](formulation::FullMatrix& cell_matrix,
          const domain::CellPtr<dim>& cell_ptr) -> void {
        formulation_ptr_->FillCellCollisionTerm(cell_matrix, cell_ptr, group);
      };
  auto boundary_bilinear_term_function =
      [&](formulation::FullMatrix& cell_matrix,
          const domain::FaceIndex face_index,
          const domain::CellPtr<dim>& cell_ptr) -> void {
    formulation_ptr_->FillBoundaryBilinearTerm(cell_matrix, cell_ptr, face_index, quadrature_point_ptr, group);
  };
  *fixed_matrix_ptr = 0;
  stamper_ptr_->StampMatrix(*fixed_matrix_ptr, streaming_term_function);
  stamper_ptr_->StampMatrix(*fixed_matrix_ptr, collision_term_function);
  stamper_ptr_->StampBoundaryMatrix(*fixed_matrix_ptr,
                                    boundary_bilinear_term_function);
}

template<int dim>
void SAAFUpdater<dim>::UpdateFissionSource(system::System &to_update,
                                           system::EnergyGroup group,
                                           quadrature::QuadraturePointIndex index) {
  auto fission_source_ptr =
      to_update.right_hand_side_ptr_->GetVariableTermPtr({group.get(), index.get()},
                                                         system::terms::VariableLinearTerms::kFissionSource);
  auto quadrature_point_ptr = quadrature_set_ptr_->GetQuadraturePoint(index);
  const auto& current_moments = to_update.current_moments->moments();
  const auto& in_group_moment = current_moments.at({group.get(), 0, 0});
  auto fission_source_function =
      [&](formulation::Vector& cell_vector,
          const domain::CellPtr<dim> &cell_ptr) -> void {
        formulation_ptr_->FillCellFissionSourceTerm(cell_vector,
                                                    cell_ptr,
                                                    quadrature_point_ptr,
                                                    group,
                                                    to_update.k_effective.value(),
                                                    in_group_moment,
                                                    current_moments);
      };
  *fission_source_ptr = 0;
  stamper_ptr_->StampVector(*fission_source_ptr, fission_source_function);
}

template<int dim>
void SAAFUpdater<dim>::UpdateScatteringSource(
    system::System &to_update,
    system::EnergyGroup group,
    quadrature::QuadraturePointIndex index) {
  auto scattering_source_ptr =
      to_update.right_hand_side_ptr_->GetVariableTermPtr({group.get(), index.get()},
                                                         system::terms::VariableLinearTerms::kScatteringSource);
  *scattering_source_ptr = 0;
  auto quadrature_point_ptr = quadrature_set_ptr_->GetQuadraturePoint(index);
  const auto& current_moments = to_update.current_moments->moments();
  const auto& in_group_moment = current_moments.at({group.get(), 0, 0});
  auto scattering_source_function =
      [&](formulation::Vector& cell_vector,
          const domain::CellPtr<dim> &cell_ptr) -> void {
    formulation_ptr_->FillCellScatteringSourceTerm(cell_vector,
                                                   cell_ptr,
                                                   quadrature_point_ptr,
                                                   group,
                                                   in_group_moment,
                                                   current_moments);
  };
  stamper_ptr_->StampVector(*scattering_source_ptr, scattering_source_function);
}

template class SAAFUpdater<1>;
template class SAAFUpdater<2>;
template class SAAFUpdater<3>;

} // namespace updater

} // namespace formulation

} // namespace bart
