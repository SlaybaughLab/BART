#ifndef BART_SRC_DOMAIN_DEFINITION_H_
#define BART_SRC_DOMAIN_DEFINITION_H_

#include <memory>

#include <deal.II/base/index_set.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/lac/constraint_matrix.h>

#include "../domain/mesh_i.h"
#include "../domain/finite_element_i.h"

namespace bart {

namespace domain {

/*! \brief Defines a domain that couples a cartesian mesh with a finite element basis.
 *
 * This class provides the framework on which the domain can be solved. The
 * purpose of the coupling is to distribute the finite element basis functions
 * and needed degrees of freedom throughout the provided mesh. Most of this is
 * done under the hood by a dealii object called the DofHandler. This class
 * is a mediator for all the required interactions between Dealii objects. It
 * maintains ownership of the mesh and finite element basis through unique pointers.
 *
 * Many portions of this class are adapted from the BARTDriver class written by
 * Weixiong Zheng.
 * 
 * \author Joshua Rehak
 * \date 2019/02
 */  

template <int dim>
class Definition {
 public:
  /*! \brief Constructor.
   * Takes ownership of injected dependencies (MeshI and FiniteElementI).
   */
  Definition(std::unique_ptr<domain::MeshI<dim>> &mesh,
             std::unique_ptr<domain::FiniteElementI<dim>> &finite_element);
  ~Definition() = default;
  
  /*! Fills triangulation with mesh defined in MeshI object
   * Creates mesh shape, sets up boundary ids and material ids. Requires that 
   * the mesh has a material mapping setup.
   */ 
  Definition<dim>& SetUpMesh();

  /*! Set up the DOF handler, to access sparsity patterns, etc */
  Definition<dim>& SetUpDOF();

  dealii::IndexSet locally_owned_dofs() const { return locally_owned_dofs_; }

  dealii::IndexSet locally_relevant_dofs() const {
    return locally_relevant_dofs_; }

  dealii::ConstraintMatrix& constraint_matrix() {
    return constraint_matrix_; }

  /*! Get total degrees of freedom */
  int total_degrees_of_freedom() const;
  
 private:
  //! Internal owned mesh object.
  std::unique_ptr<domain::MeshI<dim>> mesh_;
  
  //! Internal owned finite element object
  std::unique_ptr<domain::FiniteElementI<dim>> finite_element_;

  //! Internal distributed triangulation object
  dealii::parallel::distributed::Triangulation<dim> triangulation_;

  //! Internal DoFHandler object
  dealii::DoFHandler<dim> dof_handler_;

  //! Total degrees of freedom
  mutable int total_degrees_of_freedom_ = 0;

  //! Index of locally owned dofs owned by the current processor
  dealii::IndexSet locally_owned_dofs_;

  /*! Index of locally relevant dofs to the current processor
   * This includes indices of dofs owned by the current processor but also
   * the indices of neighboring cells (ghost cells) that may be relevant.
   */
  dealii::IndexSet locally_relevant_dofs_;

  /*! Constraint matrix */
  dealii::ConstraintMatrix constraint_matrix_;
};

} // namespace domain

} // namespace bart

#endif // BART_SRC_DOMAIN_DEFINITION_H_
