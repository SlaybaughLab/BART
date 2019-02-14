#include "domain.h"

#include <deal.II/dofs/dof_tools.h>

namespace bart {

namespace problem {

template <int dim>
Domain<dim>::Domain(std::unique_ptr<domain::MeshI<dim>> &mesh,
                    std::unique_ptr<domain::FiniteElementI<dim>> &finite_element)
    : mesh_(std::move(mesh)),                     
      finite_element_(std::move(finite_element)), 
      triangulation_(MPI_COMM_WORLD,
                     typename dealii::Triangulation<dim>::MeshSmoothing(
                         dealii::Triangulation<dim>::smoothing_on_refinement |
                         dealii::Triangulation<dim>::smoothing_on_coarsening)),
      dof_handler_(triangulation_) {}

template <int dim>
Domain<dim>& Domain<dim>::SetUpMesh() {
  // SetUp triangulation using the mesh dependency
  AssertThrow(mesh_->has_material_mapping(),
                    dealii::ExcMessage("Mesh object must have initialized material mapping"));
  mesh_->FillTriangulation(triangulation_);
  mesh_->FillBoundaryID(triangulation_);
  mesh_->FillMaterialID(triangulation_);
  return *this;
}

template <int dim>
Domain<dim>& Domain<dim>::SetUpDOF() {
  // Setup dof Handler
  dof_handler_.distribute_dofs(*(finite_element_->finite_element()));
  // Populate dof IndexSets
  locally_owned_dofs_ = dof_handler_.locally_owned_dofs();
  dealii::DoFTools::extract_locally_relevant_dofs(dof_handler_,
                                                  locally_relevant_dofs_);    
  return *this;
}

template class Domain<1>;
template class Domain<2>;

} // namespace bart

} // namespace problem
