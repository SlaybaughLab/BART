#ifndef BART_SRC_DOMAIN_DEFINITION_I_H_
#define BART_SRC_DOMAIN_DEFINITION_I_H_

#include <vector>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include "data/matrix_parameters.h"
#include "domain/domain_types.h"
#include "problem/parameter_types.h"

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
class DefinitionI {
 public:
  using CellRange = std::vector<domain::CellPtr<dim>>;

  virtual ~DefinitionI() = default;

  /*! Set up the DOF handler, to access sparsity patterns, etc */
  virtual DefinitionI<dim>& SetUpDOF() = 0;

/*! Fills triangulation with mesh defined in MeshI object
   * Creates mesh shape, sets up boundary ids and material ids. Requires that
   * the mesh has a material mapping setup.
   */
  virtual DefinitionI<dim>& SetUpMesh() = 0;

  /*! Get the parameters required to build a system matrix for this domain */
  virtual void FillMatrixParameters(
      data::MatrixParameters &to_fill,
      problem::DiscretizationType discretization) const = 0;

  /*! Get a matrix suitible for a cell matrix.
   *
   * \return a dealii FullMatrix<double> of appropriate size.
   */
  virtual dealii::FullMatrix<double> GetCellMatrix() const = 0;
  /*! Get a matrix suitible for a cell rhs.
   *
   * \return a dealii Vector<double> of appropriate size.
   */
  virtual dealii::Vector<double> GetCellVector() const = 0;

  /*! Get a range of all cells to allow iterating over them */
  virtual CellRange Cells() const = 0;

  /*! Get discretization type */
  virtual problem::DiscretizationType discretization_type() const = 0;

  /*! Get locally owned degrees of freedom */
  virtual dealii::IndexSet locally_owned_dofs() const = 0;

  /*! Get internal DOF object */
  virtual const dealii::DoFHandler<dim>& dof_handler() const = 0;

  /*! Get total degrees of freedom */
  virtual int total_degrees_of_freedom() const = 0;
};

} // namespace domain

} // namespace bart

#endif // BART_SRC_DOMAIN_DEFINITION_I_H_