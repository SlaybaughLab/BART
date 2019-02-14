#ifndef BART_SRC_PROBLEM_DOMAIN_H_
#define BART_SRC_PROBLEM_DOMAIN_H_

#include <memory>

#include <deal.II/base/index_set.h>
#include <deal.II/distributed/tria.h>

#include "../domain/mesh_i.h"
#include "../domain/finite_element_i.h"

namespace bart {

namespace problem {

/*! \brief Defines a domain that couples a cartesian mesh with a finite element basis.
 *
 * This class provides the framework on which the problem can be solved. The
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
class Domain {
 public:
  /*! \brief Constructor.
   * Takes ownership of injected dependencies (MeshI and FiniteElementI).
   */
  Domain(std::unique_ptr<domain::MeshI<dim>> &mesh,
         std::unique_ptr<domain::FiniteElementI<dim>> &finite_element);
  ~Domain() = default;
  
  /*! Fills triangulation with mesh defined in MeshI object
   * Creates mesh shape, sets up boundary ids and material ids. Requires that 
   * the mesh has a material mapping setup.
   */ 
  Domain<dim>& SetUpMesh();

  /*! Set up the DOF handler, to access sparsity patterns, etc */
  Domain<dim>& SetUpDOF();

  dealii::IndexSet locally_owned_dofs() { return locally_owned_dofs_; }
  
 private:
  //! Internal owned mesh object.
  std::unique_ptr<domain::MeshI<dim>> mesh_;
  
  //! Internal owned finite element object
  std::unique_ptr<domain::FiniteElementI<dim>> finite_element_;

  //! Internal distributed triangulation object
  dealii::parallel::distributed::Triangulation<dim> triangulation_;

  //! Internal DoFHandler object
  dealii::DoFHandler<dim> dof_handler_;

  //! Index of locally owned dofs owned by the current processor
  dealii::IndexSet locally_owned_dofs_;

  /*! Index of locally relevant dofs to the current processor
   * This includes indices of dofs owned by the current processor but also
   * the indices of neighboring cells (ghost cells) that may be relevant.
   */
  dealii::IndexSet locally_relevant_dofs_;
};

} // namespace bart

} // namespace problem

#endif // BART_SRC_PROBLEM_DOMAIN_H_
