#include "formulation/stamper.hpp"

namespace bart::formulation {

template<int dim>
Stamper<dim>::Stamper(std::shared_ptr<domain::DomainI<dim>> domain_ptr)
    : domain_ptr_(domain_ptr) {
  this->template AssertPointerNotNull(domain_ptr_.get(), "domain", "Stamper::Stamper constructor");
  this->set_description("system matrix stamper", utility::DefaultImplementation(true));
}

template<int dim>
auto Stamper<dim>::StampMatrix(system::MPISparseMatrix& to_stamp, CellMatrixStampFunction stamp_function) -> void {
  auto cell_matrix = domain_ptr_->GetCellMatrix();
  auto cells = domain_ptr_->Cells();
  std::vector<dealii::types::global_dof_index> local_dof_indices(cell_matrix.n_cols());

  for (const auto& cell : cells) {
    cell_matrix = 0;
    cell->get_dof_indices(local_dof_indices);
    stamp_function(cell_matrix, cell);
    to_stamp.add(local_dof_indices, local_dof_indices, cell_matrix);
  }
  to_stamp.compress(dealii::VectorOperation::add);
}
template<int dim>
auto Stamper<dim>::StampVector(system::MPIVector& to_stamp, CellVectorStampFunction stamp_function) -> void {
  auto cell_vector = domain_ptr_->GetCellVector();
  auto cells = domain_ptr_->Cells();
  std::vector<dealii::types::global_dof_index> local_dof_indices(cell_vector.size());

  for (const auto& cell : cells) {
    cell_vector = 0;
    cell->get_dof_indices(local_dof_indices);
    stamp_function(cell_vector, cell);
    to_stamp.add(local_dof_indices, cell_vector);
  }
  to_stamp.compress(dealii::VectorOperation::add);
}

template<int dim>
auto Stamper<dim>::StampBoundaryMatrix(system::MPISparseMatrix &to_stamp,
                                       FaceMatrixStampFunction stamp_function) -> void {
  auto cell_matrix = domain_ptr_->GetCellMatrix();
  auto cells = domain_ptr_->Cells();
  std::vector<dealii::types::global_dof_index> local_dof_indices(cell_matrix.n_cols());

  for (const auto& cell : cells) {
    if (cell->at_boundary()) {
      int faces_per_cell = dealii::GeometryInfo<dim>::faces_per_cell;
      for (int face = 0; face < faces_per_cell; ++face) {
        if (cell->face(face)->at_boundary()) {
          cell_matrix = 0;
          cell->get_dof_indices(local_dof_indices);
          stamp_function(cell_matrix, domain::FaceIndex(face), cell);
          to_stamp.add(local_dof_indices, local_dof_indices, cell_matrix);
        }
      }
    }
  }
  to_stamp.compress(dealii::VectorOperation::add);
}

template<int dim>
auto Stamper<dim>::StampBoundaryVector(system::MPIVector &to_stamp, FaceVectorStampFunction stamp_function) -> void {
  auto cell_vector = domain_ptr_->GetCellVector();
  auto cells = domain_ptr_->Cells();

  for (const auto& cell : cells) {
    if (cell->at_boundary()) {
      int faces_per_cell = dealii::GeometryInfo<dim>::faces_per_cell;
      for (int face = 0; face < faces_per_cell; ++face) {
        if (cell->face(face)->at_boundary()) {
          cell_vector = 0;
          std::vector<dealii::types::global_dof_index> local_dof_indices(
              cell_vector.size());
          cell->get_dof_indices(local_dof_indices);
          stamp_function(cell_vector, domain::FaceIndex(face), cell);
          to_stamp.add(local_dof_indices, cell_vector);
        }
      }
    }
  }
  to_stamp.compress(dealii::VectorOperation::add);
}

template class Stamper<1>;
template class Stamper<2>;
template class Stamper<3>;

} // namespace bart::formulation