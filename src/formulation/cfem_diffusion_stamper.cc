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

  auto streaming_function =
      [&](dealii::FullMatrix<double>& matrix,
          const Cell& cell_ptr) -> void
      {
        this->diffusion_ptr_->FillCellStreamingTerm(matrix,
                                                    this->diffusion_init_token_,
                                                    cell_ptr,
                                                    group);
      };

  StampMatrix(to_stamp, streaming_function);
}

template<int dim>
void CFEM_DiffusionStamper<dim>::StampCollisionTerm(MPISparseMatrix &to_stamp,
                                                    GroupNumber group) {

  auto collision_function =
      [&](dealii::FullMatrix<double>& matrix,
          const Cell& cell_ptr) -> void
      {
        this->diffusion_ptr_->FillCellCollisionTerm(matrix,
                                                    this->diffusion_init_token_,
                                                    cell_ptr,
                                                    group);
      };

  StampMatrix(to_stamp, collision_function);
}

template<int dim>
void CFEM_DiffusionStamper<dim>::StampBoundaryTerm(MPISparseMatrix &to_stamp) {

  int faces_per_cell = dealii::GeometryInfo<dim>::faces_per_cell;
  auto boundary_function =
      [&, faces_per_cell](dealii::FullMatrix<double>& matrix,
          const Cell& cell_ptr) -> void
      {
        if (cell_ptr->at_boundary()) {
          for (int face = 0; face < faces_per_cell; ++face) {
            if (cell_ptr->face(face)->at_boundary()) {
              diffusion_ptr_->FillBoundaryTerm(matrix,
                                               diffusion_init_token_,
                                               cell_ptr,
                                               face,
                                               BoundaryType::kVacuum);
            }
          }
        }
      };
  StampMatrix(to_stamp, boundary_function);
}

template <int dim>
void CFEM_DiffusionStamper<dim>::StampMatrix(
    MPISparseMatrix &to_stamp,
    std::function<void(dealii::FullMatrix<double>&, const Cell&)> function) {

  auto cell_matrix = definition_ptr_->GetCellMatrix();
  std::vector<dealii::types::global_dof_index> local_dof_indices(cell_matrix.n_cols());

  for (const auto& cell : cells_) {
    cell_matrix = 0;
    cell->get_dof_indices(local_dof_indices);
    function(cell_matrix, cell);
    to_stamp.add(local_dof_indices, local_dof_indices, cell_matrix);
  }
  to_stamp.compress(dealii::VectorOperation::add);
}


template class CFEM_DiffusionStamper<1>;
template class CFEM_DiffusionStamper<2>;
template class CFEM_DiffusionStamper<3>;

} // namespace formulation

} // namespace bart