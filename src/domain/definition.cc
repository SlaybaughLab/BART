#include "definition.h"

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/sparsity_tools.h>

namespace bart {

namespace domain {

template <int dim>
Definition<dim>::Definition(std::unique_ptr<domain::MeshI<dim>> &mesh,
                    std::unique_ptr<domain::FiniteElementI<dim>> &finite_element)
    : mesh_(std::move(mesh)),                     
      finite_element_(std::move(finite_element)), 
      triangulation_(MPI_COMM_WORLD,
                     typename dealii::Triangulation<dim>::MeshSmoothing(
                         dealii::Triangulation<dim>::smoothing_on_refinement |
                         dealii::Triangulation<dim>::smoothing_on_coarsening)),
      dof_handler_(triangulation_) {}

template <int dim>
Definition<dim>& Definition<dim>::SetUpMesh() {
  // SetUp triangulation using the mesh dependency
  AssertThrow(mesh_->has_material_mapping(),
                    dealii::ExcMessage("Mesh object must have initialized material mapping"));
  mesh_->FillTriangulation(triangulation_);
  mesh_->FillBoundaryID(triangulation_);
  mesh_->FillMaterialID(triangulation_);
  return *this;
}

template <int dim>
Definition<dim>& Definition<dim>::SetUpDOF() {
  // Setup dof Handler
  dof_handler_.distribute_dofs(*(finite_element_->finite_element()));
  // Populate dof IndexSets
  locally_owned_dofs_ = dof_handler_.locally_owned_dofs();
  dealii::DoFTools::extract_locally_relevant_dofs(dof_handler_,
                                                  locally_relevant_dofs_);
  // Create constraint matrix
  constraint_matrix_.clear();
  constraint_matrix_.reinit(locally_relevant_dofs_);
  dealii::DoFTools::make_hanging_node_constraints(dof_handler_,
                                                  constraint_matrix_);
  constraint_matrix_.close();
  return *this;
}

template <int dim>
void Definition<dim>::FillMatrixParameters(
    data::MatrixParameters &to_fill,
    problem::DiscretizationType discretization) const {

  to_fill.row_degrees_of_freedom = locally_owned_dofs_;
  to_fill.column_degrees_of_freedom = locally_owned_dofs_;

  dealii::DynamicSparsityPattern dsp(locally_relevant_dofs_);

  if (discretization ==  problem::DiscretizationType::kDiscontinuousFEM) {
    dealii::DoFTools::make_flux_sparsity_pattern(dof_handler_, dsp,
                                                 constraint_matrix_, false);
  } else {
    dealii::DoFTools::make_sparsity_pattern(dof_handler_, dsp,
                                            constraint_matrix_, false);
  }

  auto locally_owned_per_proc =
      dof_handler_.n_locally_owned_dofs_per_processor();
  
  dealii::SparsityTools::distribute_sparsity_pattern(dsp,
                                                     locally_owned_per_proc,
                                                     MPI_COMM_WORLD,
                                                     locally_relevant_dofs_);
  to_fill.dsp = dsp;
}

template <int dim>
int Definition<dim>::total_degrees_of_freedom() const {
  if (total_degrees_of_freedom_ == 0)
    total_degrees_of_freedom_ = dof_handler_.n_dofs();
  return total_degrees_of_freedom_;
}

template class Definition<1>;
template class Definition<2>;

} // namespace bart

} // namespace domain
