#include "formulation/cfem_diffusion_stamper.h"
#include "cfem_diffusion_stamper.h"

namespace bart {

namespace formulation {

template<int dim>
CFEM_DiffusionStamper<dim>::CFEM_DiffusionStamper(
    std::unique_ptr<formulation::scalar::DiffusionI<dim>> diffusion_ptr,
    std::shared_ptr<domain::DefinitionI<dim>> definition_ptr,
    const std::unordered_set<Boundary> reflective_boundaries)
      : diffusion_ptr_(std::move(diffusion_ptr)),
        definition_ptr_(definition_ptr),
        reflective_boundaries_(reflective_boundaries) {

  cells_ = definition_ptr_->Cells();
  diffusion_ptr_->Precalculate(cells_[0]);
}

template<int dim>
CFEM_DiffusionStamper<dim>::CFEM_DiffusionStamper(
    std::unique_ptr<formulation::scalar::DiffusionI<dim>> diffusion_ptr,
    std::shared_ptr<domain::DefinitionI<dim>> definition_ptr,
    const std::map<Boundary, bool> reflective_boundary_map)
    : CFEM_DiffusionStamper<dim>(std::move(diffusion_ptr),
                                 definition_ptr) {
  for (auto const pair : reflective_boundary_map) {
    auto [boundary, is_reflective] = pair;
    if (is_reflective)
      reflective_boundaries_.insert(boundary);
  }
}


template<int dim>
void CFEM_DiffusionStamper<dim>::StampStreamingTerm(MPISparseMatrix &to_stamp,
                                                    GroupNumber group) {

  auto streaming_function =
      [&](dealii::FullMatrix<double>& matrix,
          const domain::CellPtr<dim>& cell_ptr) -> void
      {
        this->diffusion_ptr_->FillCellStreamingTerm(matrix,
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
          const domain::CellPtr<dim>& cell_ptr) -> void
      {
        this->diffusion_ptr_->FillCellCollisionTerm(matrix,
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
          const domain::CellPtr<dim>& cell_ptr) -> void
      {
        if (cell_ptr->at_boundary()) {
          for (int face = 0; face < faces_per_cell; ++face) {
            if (cell_ptr->face(face)->at_boundary()) {

              Boundary boundary = static_cast<Boundary>(
                  cell_ptr->face(face)->boundary_id());
              BoundaryType boundary_type = BoundaryType::kVacuum;

              if (reflective_boundaries_.count(boundary) == 1)
                boundary_type = BoundaryType::kReflective;

              diffusion_ptr_->FillBoundaryTerm(matrix,
                                               cell_ptr,
                                               face,
                                               boundary_type);
            }
          }
        }
      };
  StampMatrix(to_stamp, boundary_function);
}

template<int dim>
void CFEM_DiffusionStamper<dim>::StampFixedSource(MPIVector &to_stamp,
                                                  const GroupNumber group) {
  auto fixed_source_function =
      [&](dealii::Vector<double> &vector,
          const domain::CellPtr<dim>& cell_ptr) -> void {
    this->diffusion_ptr_->FillCellFixedSource(vector, cell_ptr, group);
  };
  StampVector(to_stamp, fixed_source_function);
}

template<int dim>
void CFEM_DiffusionStamper<dim>::StampFissionSource(
    MPIVector& to_stamp,
    const GroupNumber group,
    const double k_effective,
    const system::moments::MomentVector& in_group_moment,
    const system::moments::MomentsMap& group_moments) {

  auto fission_source_function =
      [&](dealii::Vector<double> &vector,
          const domain::CellPtr<dim>& cell_ptr) -> void {
    this->diffusion_ptr_->FillCellFissionSource(vector,
                                                cell_ptr,
                                                group,
                                                k_effective,
                                                in_group_moment,
                                                group_moments);
  };
  StampVector(to_stamp, fission_source_function);
}

template<int dim>
void CFEM_DiffusionStamper<dim>::StampScatteringSource(
    MPIVector &to_stamp,
    const GroupNumber group,
    const system::moments::MomentsMap &group_moments) {
  auto scattering_source_function =
      [&](dealii::Vector<double> &vector,
          const domain::CellPtr<dim>& cell_ptr) -> void {
        this->diffusion_ptr_->FillCellScatteringSource(vector,
                                                       cell_ptr,
                                                       group,
                                                       group_moments);
      };
  StampVector(to_stamp, scattering_source_function);
}

template <int dim>
void CFEM_DiffusionStamper<dim>::StampMatrix(
    MPISparseMatrix &to_stamp,
    std::function<void(dealii::FullMatrix<double>&,
        const domain::CellPtr<dim>&)> function) {

  auto cell_matrix = definition_ptr_->GetCellMatrix();
  std::vector<dealii::types::global_dof_index> local_dof_indices(
      cell_matrix.n_cols());

  for (const auto& cell : cells_) {
    cell_matrix = 0;
    cell->get_dof_indices(local_dof_indices);
    function(cell_matrix, cell);
    to_stamp.add(local_dof_indices, local_dof_indices, cell_matrix);
  }
  to_stamp.compress(dealii::VectorOperation::add);
}
template<int dim>
void CFEM_DiffusionStamper<dim>::StampVector(
    MPIVector &to_stamp,
    std::function<void(dealii::Vector<double> &,
                       const domain::CellPtr<dim> &)> function) {

  auto cell_vector = definition_ptr_->GetCellVector();
  std::vector<dealii::types::global_dof_index> local_dof_indices(
      cell_vector.size());

  for (const auto& cell : cells_) {
    cell_vector = 0;
    cell->get_dof_indices(local_dof_indices);
    function(cell_vector, cell);
    to_stamp.add(local_dof_indices, cell_vector);
  }
  to_stamp.compress(dealii::VectorOperation::add);
}


template class CFEM_DiffusionStamper<1>;
template class CFEM_DiffusionStamper<2>;
template class CFEM_DiffusionStamper<3>;

} // namespace formulation

} // namespace bart