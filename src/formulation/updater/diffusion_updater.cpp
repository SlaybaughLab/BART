#include "formulation/updater/diffusion_updater.hpp"

namespace bart {

namespace formulation {

namespace updater {

template<int dim>
DiffusionUpdater<dim>::DiffusionUpdater(
    std::unique_ptr<DiffusionFormulationType> formulation_ptr,
    std::unique_ptr<StamperType> stamper_ptr,
    std::unordered_set<problem::Boundary> reflective_boundaries)
    : formulation_ptr_(std::move(formulation_ptr)),
      stamper_ptr_(std::shared_ptr<StamperType>(std::move(stamper_ptr))),
      reflective_boundaries_(reflective_boundaries) {
  AssertThrow(formulation_ptr_ != nullptr,
              dealii::ExcMessage("Error in constructor of DiffusionUpdater, "
                                 "formulation pointer passed is null"))
  AssertThrow(stamper_ptr_ != nullptr,
              dealii::ExcMessage("Error in constructor of DiffusionUpdater, "
                                 "stamper pointer passed is null"))
  this->set_description("diffusion formulation updater",
                        utility::DefaultImplementation(true));

  if (!reflective_boundaries_.empty()) {
    this->set_description(this->description() + " (reflective BCs)",
                          utility::DefaultImplementation(true));
  }
}

template<int dim>
void DiffusionUpdater<dim>::UpdateFixedTerms(
    system::System& to_update,
    system::EnergyGroup energy_group,
    quadrature::QuadraturePointIndex /*index*/) {
  using CellPtr = domain::CellPtr<dim>;
  int group = energy_group.get();
  auto fixed_matrix_ptr =
     to_update.left_hand_side_ptr_->GetFixedTermPtr({group, 0});
  auto fixed_vector_ptr =
      to_update.right_hand_side_ptr_->GetFixedTermPtr({group, 0});
  auto streaming_term_function = [&](formulation::FullMatrix& cell_matrix,
                                     const CellPtr& cell_ptr) -> void {
        formulation_ptr_->FillCellStreamingTerm(cell_matrix, cell_ptr, group);
      };
  auto collision_term_function = [&](formulation::FullMatrix& cell_matrix,
                                     const CellPtr& cell_ptr) -> void {
        formulation_ptr_->FillCellCollisionTerm(cell_matrix, cell_ptr, group);
      };
  auto fixed_term_function = [&](formulation::Vector& cell_vector,
                                 const CellPtr& cell_ptr) -> void {
    formulation_ptr_->FillCellFixedSource(cell_vector, cell_ptr, group);
  };
  auto boundary_function = [&](formulation::FullMatrix& cell_matrix,
                               const domain::FaceIndex face_index,
                               const CellPtr& cell_ptr) -> void {
    using DiffusionBoundaryType = typename formulation::scalar::DiffusionI<dim>::BoundaryType;
    problem::Boundary boundary = static_cast<problem::Boundary>(
        cell_ptr->face(face_index.get())->boundary_id());
    DiffusionBoundaryType boundary_type = DiffusionBoundaryType::kVacuum;

    if (reflective_boundaries_.count(boundary) == 1)
      boundary_type = DiffusionBoundaryType::kReflective;
    formulation_ptr_->FillBoundaryTerm(cell_matrix, cell_ptr,
                                       face_index.get(), boundary_type);
  };
  *fixed_matrix_ptr = 0;
  *fixed_vector_ptr = 0;
  stamper_ptr_->StampMatrix(*fixed_matrix_ptr, streaming_term_function);
  stamper_ptr_->StampMatrix(*fixed_matrix_ptr, collision_term_function);
  stamper_ptr_->StampBoundaryMatrix(*fixed_matrix_ptr, boundary_function);
  stamper_ptr_->StampVector(*fixed_vector_ptr, fixed_term_function);
}
template<int dim>
void DiffusionUpdater<dim>::UpdateScatteringSource(
    system::System &to_update,
    system::EnergyGroup energy_group,
    quadrature::QuadraturePointIndex /*index*/) {
  int group = energy_group.get();
  auto scattering_source_ptr =
      to_update.right_hand_side_ptr_->GetVariableTermPtr({group, 0},
                                                          system::terms::VariableLinearTerms::kScatteringSource);
  *scattering_source_ptr = 0;
  const auto& current_moments = to_update.current_moments->moments();
  auto scattering_source_function =
      [&](formulation::Vector& cell_vector,
          const domain::CellPtr<dim> &cell_ptr) -> void {
        formulation_ptr_->FillCellScatteringSource(cell_vector,
                                                   cell_ptr,
                                                   group,
                                                   current_moments);
      };
  stamper_ptr_->StampVector(*scattering_source_ptr, scattering_source_function);
}
template<int dim>
void DiffusionUpdater<dim>::UpdateFissionSource(
    system::System &to_update,
    system::EnergyGroup energy_group,
    quadrature::QuadraturePointIndex /*index*/) {
  int group = energy_group.get();
  auto scattering_source_ptr =
      to_update.right_hand_side_ptr_->GetVariableTermPtr({group, 0},
                                                         system::terms::VariableLinearTerms::kFissionSource);
  *scattering_source_ptr = 0;
  const auto& current_moments = to_update.current_moments->moments();
  const auto& in_group_moment = current_moments.at({group, 0, 0});
  auto fission_source_function =
      [&](formulation::Vector& cell_vector,
          const domain::CellPtr<dim> &cell_ptr) -> void {
        formulation_ptr_->FillCellFissionSource(cell_vector,
                                                cell_ptr,
                                                group,
                                                to_update.k_effective.value(),
                                                in_group_moment,
                                                current_moments);
      };
  stamper_ptr_->StampVector(*scattering_source_ptr, fission_source_function);
}
template<int dim>
void DiffusionUpdater<dim>::UpdateFixedSource(
    system::System &to_update,
    system::EnergyGroup energy_group,
    quadrature::QuadraturePointIndex /*index*/) {
  int group = energy_group.get();
  auto fixed_source_ptr =
      to_update.right_hand_side_ptr_->GetFixedTermPtr({group, 0});
  *fixed_source_ptr = 0;
  auto fixed_source_function =
      [&](formulation::Vector& cell_vector,
          const domain::CellPtr<dim> &cell_ptr) -> void {
        formulation_ptr_->FillCellFixedSource(cell_vector,
                                              cell_ptr, group);
      };
  stamper_ptr_->StampVector(*fixed_source_ptr, fixed_source_function);
}

template class DiffusionUpdater<1>;
template class DiffusionUpdater<2>;
template class DiffusionUpdater<3>;

} // namespace updater

} // namespace formulation

} // namespace bart
