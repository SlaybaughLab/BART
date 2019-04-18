#include "formulation/cfem_diffusion_stamper.h"
#include "cfem_diffusion_stamper.h"

namespace bart {

namespace formulation {

template<int dim>
CFEM_DiffusionStamper<dim>::CFEM_DiffusionStamper(
    std::unique_ptr<formulation::scalar::CFEM_DiffusionI<dim>> diffusion_ptr,
    std::unique_ptr<domain::DefinitionI<dim>> definition_ptr)
      : diffusion_ptr_(std::move(diffusion_ptr)),
        definition_ptr_(std::move(definition_ptr)) {

  cells_ = definition_ptr_->Cells();
  diffusion_init_token_ = diffusion_ptr_->Precalculate(cells_[0]);
}



template<int dim>
void CFEM_DiffusionStamper<dim>::StampStreamingTerm(MPISparseMatrix &to_stamp,
                                                    GroupNumber group) {

  auto cell_matrix = definition_ptr_->GetCellMatrix();
  int material_id = 0;
  std::vector<dealii::types::global_dof_index> local_dof_indices(cell_matrix.n_cols());

  for (const auto& cell : cells_) {
    material_id = cell->material_id();
    cell->get_dof_indices(local_dof_indices);
    diffusion_ptr_->FillCellStreamingTerm(cell_matrix, diffusion_init_token_, cell,
                                          material_id, group);
    to_stamp.add(local_dof_indices, local_dof_indices, cell_matrix);
  }
  to_stamp.compress(dealii::VectorOperation::add);
}

template<int dim>
void CFEM_DiffusionStamper<dim>::StampCollisionTerm(MPISparseMatrix &to_stamp,
                                                    GroupNumber group) {

  auto cell_matrix = definition_ptr_->GetCellMatrix();
  int material_id = 0;
  std::vector<dealii::types::global_dof_index> local_dof_indices(cell_matrix.n_cols());

  for (const auto& cell : cells_) {
    material_id = cell->material_id();
    cell->get_dof_indices(local_dof_indices);
    diffusion_ptr_->FillCellCollisionTerm(cell_matrix, diffusion_init_token_, cell,
                                          material_id, group);
    to_stamp.add(local_dof_indices, local_dof_indices, cell_matrix);
  }
  to_stamp.compress(dealii::VectorOperation::add);
}

template class CFEM_DiffusionStamper<1>;
template class CFEM_DiffusionStamper<2>;
template class CFEM_DiffusionStamper<3>;

} // namespace formulation

} // namespace bart