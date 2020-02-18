#include "definition.h"

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

namespace bart {

namespace domain {

template <int dim>
Definition<dim>::Definition(
    std::unique_ptr<domain::mesh::MeshI<dim>> mesh,
    std::shared_ptr<domain::finite_element::FiniteElementI<dim>> finite_element,
    problem::DiscretizationType discretization)
    : mesh_(std::move(mesh)),                     
      finite_element_(finite_element),
      triangulation_(MPI_COMM_WORLD,
                     typename dealii::Triangulation<dim>::MeshSmoothing(
                         dealii::Triangulation<dim>::smoothing_on_refinement |
                         dealii::Triangulation<dim>::smoothing_on_coarsening)),
      dof_handler_(triangulation_),
      discretization_type_(discretization) {}

template <>
Definition<1>::Definition(
    std::unique_ptr<domain::mesh::MeshI<1>> mesh,
    std::shared_ptr<domain::finite_element::FiniteElementI<1>> finite_element,
    problem::DiscretizationType discretization)
    : mesh_(std::move(mesh)),
      finite_element_(finite_element),
      triangulation_(typename dealii::Triangulation<1>::MeshSmoothing(
                         dealii::Triangulation<1>::smoothing_on_refinement |
                             dealii::Triangulation<1>::smoothing_on_coarsening)),
      dof_handler_(triangulation_),
      discretization_type_(discretization) {}

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

  for (auto cell = dof_handler_.begin_active();
       cell != dof_handler_.end(); ++cell) {
    if (cell->is_locally_owned())
      local_cells_.push_back(cell);
  }

  // Set up dynamic sparsity pattern
  dynamic_sparsity_pattern_.reinit(locally_relevant_dofs_.size(),
                                   locally_relevant_dofs_.size(),
                                   locally_relevant_dofs_);

  if (discretization_type_ ==  problem::DiscretizationType::kDiscontinuousFEM) {
    dealii::DoFTools::make_flux_sparsity_pattern(dof_handler_,
                                                 dynamic_sparsity_pattern_,
                                                 constraint_matrix_, false);
  } else {
    dealii::DoFTools::make_sparsity_pattern(dof_handler_,
                                            dynamic_sparsity_pattern_,
                                            constraint_matrix_, false);
  }

  dealii::SparsityTools::distribute_sparsity_pattern(
      dynamic_sparsity_pattern_,
      dof_handler_.n_locally_owned_dofs_per_processor(),
      MPI_COMM_WORLD, locally_relevant_dofs_);

  constraint_matrix_.condense(dynamic_sparsity_pattern_);

  return *this;
}

template <>
Definition<1>& Definition<1>::SetUpDOF() {
  auto n_mpi_processes = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  auto this_process = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  dealii::GridTools::partition_triangulation(n_mpi_processes, triangulation_);
  dof_handler_.distribute_dofs(*(finite_element_)->finite_element());
  dealii::DoFRenumbering::subdomain_wise(dof_handler_);

  for (auto cell = dof_handler_.begin_active();
       cell != dof_handler_.end(); ++cell) {
    if (cell->is_locally_owned())
      local_cells_.push_back(cell);
  }

  auto locally_owned_dofs_vector =
      dealii::DoFTools::locally_owned_dofs_per_subdomain(dof_handler_);
  locally_owned_dofs_ = locally_owned_dofs_vector.at(this_process);

  constraint_matrix_.clear();
  dealii::DoFTools::make_hanging_node_constraints(dof_handler_,
                                                  constraint_matrix_);
  constraint_matrix_.close();

  dynamic_sparsity_pattern_.reinit(dof_handler_.n_dofs(), dof_handler_.n_dofs());

  if (discretization_type_ ==  problem::DiscretizationType::kDiscontinuousFEM) {
    dealii::DoFTools::make_flux_sparsity_pattern(dof_handler_, dynamic_sparsity_pattern_,
                                                 constraint_matrix_, false);
  } else {
    dealii::DoFTools::make_sparsity_pattern(dof_handler_, dynamic_sparsity_pattern_,
                                            constraint_matrix_, false);
  }

  return *this;
}

template<int dim>
std::shared_ptr<system::MPISparseMatrix> Definition<dim>::MakeSystemMatrix() const {
  auto system_matrix_ptr = std::make_shared<system::MPISparseMatrix>();
  system_matrix_ptr->reinit(locally_owned_dofs_,
      locally_owned_dofs_,
      dynamic_sparsity_pattern_,
      MPI_COMM_WORLD);
  return system_matrix_ptr;
}

template <int dim>
void Definition<dim>::FillMatrixParameters(
    data::MatrixParameters &to_fill,
    problem::DiscretizationType discretization_type) const {

  AssertThrow(dof_handler_.has_active_dofs(),
              dealii::ExcMessage("SetUpDOF must be called before MatrixParameters are generated"));              
  
  to_fill.rows = locally_owned_dofs_;
  to_fill.columns = locally_owned_dofs_;

  dealii::DynamicSparsityPattern& dsp = to_fill.sparsity_pattern;
  dsp.reinit(locally_relevant_dofs_.size(),
             locally_relevant_dofs_.size(),
             locally_relevant_dofs_);

  if (discretization_type ==  problem::DiscretizationType::kDiscontinuousFEM) {
    dealii::DoFTools::make_flux_sparsity_pattern(dof_handler_, dsp,
                                                 constraint_matrix_, false);
  } else {
    dealii::DoFTools::make_sparsity_pattern(dof_handler_, dsp,
                                            constraint_matrix_, false);
  }

  dealii::SparsityTools::distribute_sparsity_pattern(
      dsp, dof_handler_.n_locally_owned_dofs_per_processor(),
      MPI_COMM_WORLD, locally_relevant_dofs_);

  constraint_matrix_.condense(dsp);
}

template <>
void Definition<1>::FillMatrixParameters(
    data::MatrixParameters &to_fill,
    problem::DiscretizationType discretization_type) const {

  AssertThrow(dof_handler_.has_active_dofs(),
              dealii::ExcMessage("SetUpDOF must be called before MatrixParameters are generated"));

  to_fill.rows = locally_owned_dofs_;
  to_fill.columns = locally_owned_dofs_;

  dealii::DynamicSparsityPattern& dsp = to_fill.sparsity_pattern;
  dsp.reinit(dof_handler_.n_dofs(), dof_handler_.n_dofs());

  if (discretization_type ==  problem::DiscretizationType::kDiscontinuousFEM) {
    dealii::DoFTools::make_flux_sparsity_pattern(dof_handler_, dsp,
                                                 constraint_matrix_, false);
  } else {
    dealii::DoFTools::make_sparsity_pattern(dof_handler_, dsp,
                                            constraint_matrix_, false);
  }
}

template <int dim>
int Definition<dim>::total_degrees_of_freedom() const {
  if (total_degrees_of_freedom_ == 0)
    total_degrees_of_freedom_ = dof_handler_.n_dofs();
  return total_degrees_of_freedom_;
}


template class Definition<1>;
template class Definition<2>;
template class Definition<3>;

} // namespace bart

} // namespace domain
