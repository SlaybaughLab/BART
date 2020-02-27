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

template class SAAFUpdater<1>;
template class SAAFUpdater<2>;
template class SAAFUpdater<3>;

} // namespace updater

} // namespace formulation

} // namespace bart
