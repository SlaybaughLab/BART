#include "formulation/assemble_scalar.h"

#include "formulation/types.h"

namespace bart {

namespace formulation {

template <int dim>
AssembleScalar<dim>::AssembleScalar(
    std::unique_ptr<equation::TransportScalar<dim>> equation,
    std::unique_ptr<domain::Definition<dim>> domain,
    std::shared_ptr<data::SystemScalarFluxes> scalar_fluxes,
    std::shared_ptr<data::ScalarSystemMatrices> system_matrices,
    std::shared_ptr<data::ScalarRightHandSideVectors> right_hand_side,
    std::unique_ptr<data::ScalarRightHandSideVectors> fixed_right_hand_side,
    std::map<problem::Boundary, bool> reflective_boundary_map)
    : equation_(std::move(equation)),
      domain_(std::move(domain)),
      scalar_fluxes_(scalar_fluxes),
      system_matrices_(system_matrices),
      right_hand_side_(right_hand_side),
      fixed_right_hand_side_(std::move(fixed_right_hand_side)),
      reflective_boundary_map_(reflective_boundary_map) {
  equation_->Precalculate(domain_->Cells()[0]);

}

template<int dim>
void AssembleScalar<dim>::AssembleBilinearTerms(GroupNumber group) {
  CellMatrix cell_matrix = domain_->GetCellMatrix();
  std::vector<dealii::types::global_dof_index> local_indices;
  int faces_per_cell = dealii::GeometryInfo<dim>::faces_per_cell;


  (*system_matrices_)[group] = 0.0;

  for (const auto &cell : domain_->Cells()) {
    cell_matrix = 0.0;
    equation_->FillCellBilinearTerm(cell_matrix, cell, group);

    for (int face_number = 0; face_number < faces_per_cell; ++face_number) {

      if (cell->at_boundary(face_number)) {
        auto boundary_type = BoundaryType::kVacuum;
        auto boundary_name = static_cast<problem::Boundary>(
            cell->face(face_number)->boundary_id());
        bool is_reflective = reflective_boundary_map_[boundary_name];

        if (is_reflective)
          boundary_type = BoundaryType::kReflective;

        equation_->FillBoundaryBilinearTerm(cell_matrix, cell, group,
                                                face_number, boundary_type);
      }
    }

    cell->get_dof_indices(local_indices);
    (*system_matrices_)[group].add(local_indices, local_indices, cell_matrix);
  }
  (*system_matrices_)[group].compress(dealii::VectorOperation::add);
}

template<int dim>
void AssembleScalar<dim>::AssembleFixedSourceLinearTerm(GroupNumber group) {
  CellVector cell_vector = domain_->GetCellVector();
  std::vector<dealii::types::global_dof_index> local_indices;

  (*fixed_right_hand_side_)[group] = 0.0;

  for (const auto &cell : domain_->Cells()) {
    cell_vector = 0.0;
    cell->get_dof_indices(local_indices);
    equation_->FillCellFixedSourceLinearTerm(cell_vector, cell, group);
    (*fixed_right_hand_side_)[group].add(local_indices, cell_vector);
  }

  (*fixed_right_hand_side_)[group].compress(dealii::VectorOperation::add);
}

template class AssembleScalar<1>;
template class AssembleScalar<2>;
template class AssembleScalar<3>;

} // namespace formulation

} // namespace bart

