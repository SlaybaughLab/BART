#include "formulation/cfem_saaf_stamper.h"

namespace bart {

namespace formulation {

template<int dim>
CFEM_SAAF_Stamper<dim>::CFEM_SAAF_Stamper(
    std::unique_ptr<SAAFFormulationType> saaf_ptr,
    std::shared_ptr<DomainDefinitionType> defintion_ptr)
    : formulation_ptr_(std::move(saaf_ptr)),
      definition_ptr_(defintion_ptr) {

  cells_ = definition_ptr_->Cells();
  saaf_initialization_token_ = formulation_ptr_->Initialize(cells_.at(0));
}


template<int dim>
void CFEM_SAAF_Stamper<dim>::StampCollisionTerm(
    system::MPISparseMatrix &to_stamp,
    const system::EnergyGroup group_number) {
  auto collision_term_function =
      [&](formulation::FullMatrix& cell_matrix,
          const formulation::CellPtr<dim>& cell_ptr) -> void {
        formulation_ptr_->FillCellCollisionTerm(cell_matrix,
                                                saaf_initialization_token_,
                                                cell_ptr,
                                                group_number);
  };

  StampMatrix(to_stamp, collision_term_function);
}

template<int dim>
void CFEM_SAAF_Stamper<dim>::StampMatrix(
    system::MPISparseMatrix &to_stamp,
    std::function<void(formulation::FullMatrix&,
                       const formulation::CellPtr<dim> &)> stamping_function) {
  auto cell_matrix = definition_ptr_->GetCellMatrix();
  std::vector<dealii::types::global_dof_index>
      local_dof_indices(cell_matrix.n_cols());

  for (const auto& cell : cells_) {
    cell_matrix = 0;
    cell->get_dof_indices(local_dof_indices);
    stamping_function(cell_matrix, cell);
    to_stamp.add(local_dof_indices, local_dof_indices, cell_matrix);
  }
  to_stamp.compress(dealii::VectorOperation::add);
}


template class CFEM_SAAF_Stamper<1>;
template class CFEM_SAAF_Stamper<2>;
template class CFEM_SAAF_Stamper<3>;

} // namespace formulation

} // namespace bart
