#ifndef BART_SRC_DOMAIN_DOMAIN_I_HPP_
#define BART_SRC_DOMAIN_DOMAIN_I_HPP_

#include <vector>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include "domain/domain_types.h"
#include "problem/parameter_types.hpp"
#include "system/system_types.h"
#include "utility/has_description.h"

namespace bart::domain {

/*! \brief Defines a domain that couples a cartesian mesh with a finite element basis.
 *
 * This class provides the framework on which the domain can be solved. The
 * purpose of the coupling is to distribute the finite element basis functions
 * and needed degrees of freedom throughout the provided mesh. Most of this is
 * done under the hood by a dealii object called the DofHandler. This class
 * is a mediator for all the required interactions between Dealii objects. It
 * maintains ownership of the mesh and finite element basis through unique pointers.
 *
 * The class has two major functions: (1) it acts as a wrapper for the necessary
 * dealii classes to couple a finite element basis to a mesh, and (2) provides
 * matrices based on this coupling. The functions SetUpDOF() and SetUpMesh() are
 * required for initializing the underlying objects.
 *
 * Many portions of this class are adapted from the BARTDriver class written by
 * Weixiong Zheng.
 *
 * \author Joshua Rehak
 * \date 2019/02
 */
template <int dim>
class DomainI : public utility::HasDescription {
 public:
  using CellRange = std::vector<domain::CellPtr<dim>>;

  virtual ~DomainI() = default;

  /*! Set up the DOF handler, to access sparsity patterns, etc */
  virtual auto SetUpDOF() -> DomainI<dim>& = 0;

  /*! \brief Fills triangulation with mesh defined in MeshI object
   * Creates mesh shape, sets up boundary ids and material ids. Requires that
   * the mesh has a material mapping setup.
   */
  virtual auto SetUpMesh() -> DomainI<dim>& = 0;

  /*! \brief Fills triangulation with mesh defined in MeshI object.
   *
   * Creates mesh shape, sets up boundary ids and material ids. Requires that
   * the mesh has a material mapping setup. Also performs global refinements.
   *
   * @param global_refinements number of global refinements to perform
   */
  virtual auto SetUpMesh(const int global_refinements) -> DomainI<dim>& = 0;

  /*! Get a matrix suitible for a cell matrix.
   *
   * \return a dealii FullMatrix<double> of appropriate size.
   */
  virtual auto GetCellMatrix() const -> dealii::FullMatrix<double> = 0;
  /*! Get a matrix suitible for a cell rhs.
   *
   * \return a dealii Vector<double> of appropriate size.
   */
  virtual auto GetCellVector() const -> dealii::Vector<double> = 0;

  /*! Get an MPI matrix suitable for the system */
  virtual auto MakeSystemMatrix() const -> std::shared_ptr<bart::system::MPISparseMatrix> = 0;

  /*! Get an MPI vector suitable for the system */
  virtual auto MakeSystemVector() const -> std::shared_ptr<bart::system::MPIVector> = 0;

  /*! Get a range of all cells to allow iterating over them */
  virtual auto Cells() const -> CellRange = 0;

  /*! Get discretization type */
  virtual auto discretization_type() const -> problem::DiscretizationType = 0;

  /*! Get locally owned degrees of freedom */
  virtual auto locally_owned_dofs() const -> dealii::IndexSet = 0;

  /*! Get internal DOF object */
  virtual auto dof_handler() const -> const dealii::DoFHandler<dim>& = 0;

  /*! Get total degrees of freedom */
  virtual auto total_degrees_of_freedom() const -> int = 0;
};

} // namespace bart::domain

#endif // BART_SRC_DOMAIN_DOMAIN_I_HPP_