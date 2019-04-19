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

  auto streaming_function = [&](dealii::FullMatrix<double>& matrix,
                                   const Cell& cell_ptr,
                                   const int material_id) -> void
  { this->diffusion_ptr_->FillCellStreamingTerm(matrix,
                                                this->diffusion_init_token_,
                                                cell_ptr,
                                                material_id,
                                                group);};
  StampMatrix(to_stamp, streaming_function);
}

template<int dim>
void CFEM_DiffusionStamper<dim>::StampCollisionTerm(MPISparseMatrix &to_stamp,
                                                    GroupNumber group) {

  auto collision_function = [&](dealii::FullMatrix<double>& matrix,
                                const Cell& cell_ptr,
                                const int material_id) -> void
  { this->diffusion_ptr_->FillCellCollisionTerm(matrix,
                                                this->diffusion_init_token_,
                                                cell_ptr,
                                                material_id,
                                                group);};
  StampMatrix(to_stamp, collision_function);
}

template <int dim>
void CFEM_DiffusionStamper<dim>::StampMatrix(
    MPISparseMatrix &to_stamp,
    std::function<void(dealii::FullMatrix<double>&, const Cell&, const int)> function) {

  auto cell_matrix = definition_ptr_->GetCellMatrix();
  int material_id = 0;
  std::vector<dealii::types::global_dof_index> local_dof_indices(cell_matrix.n_cols());

  for (const auto& cell : cells_) {
    material_id = cell->material_id();
    cell->get_dof_indices(local_dof_indices);
    function(cell_matrix, cell, material_id);
    to_stamp.add(local_dof_indices, local_dof_indices, cell_matrix);
  }
  to_stamp.compress(dealii::VectorOperation::add);
}

template class CFEM_DiffusionStamper<1>;
template class CFEM_DiffusionStamper<2>;
template class CFEM_DiffusionStamper<3>;

} // namespace formulation

} // namespace bart