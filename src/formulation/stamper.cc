#include "formulation/stamper.h"

namespace bart {

namespace formulation {

template<int dim>
Stamper<dim>::Stamper(std::shared_ptr<domain::DefinitionI<dim>> domain_ptr)
    : domain_ptr_(domain_ptr) {
  AssertThrow(domain_ptr_ != nullptr,
      dealii::ExcMessage("Error in constructor of formulation::Stamper, "
                         "provided domain::DefinitionI pointer is null"))
}
template<int dim>
void Stamper<dim>::StampMatrix(
    system::MPISparseMatrix& to_stamp,
    std::function<void(formulation::FullMatrix&,
                       const domain::CellPtr<dim> &)> stamp_function) {
  auto cell_matrix = domain_ptr_->GetCellMatrix();
  auto cells = domain_ptr_->Cells();
  std::vector<dealii::types::global_dof_index> local_dof_indices(
      cell_matrix.n_cols());

  for (const auto& cell : cells) {
    cell_matrix = 0;
    cell->get_dof_indices(local_dof_indices);
    stamp_function(cell_matrix, cell);
    to_stamp.add(local_dof_indices, local_dof_indices, cell_matrix);
  }
  to_stamp.compress(dealii::VectorOperation::add);
}

template class Stamper<1>;
template class Stamper<2>;
template class Stamper<3>;

} // namespace formulation

} // namespace bart
