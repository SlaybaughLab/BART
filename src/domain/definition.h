#ifndef BART_SRC_DOMAIN_DEFINITION_H_
#define BART_SRC_DOMAIN_DEFINITION_H_

#include <memory>
#include <vector>

#include <deal.II/base/index_set.h>
#include <deal.II/base/iterator_range.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/lac/constraint_matrix.h>

#include "data/matrix_parameters.h"
#include "domain/definition_i.h"
#include "domain/finite_element/finite_element_i.h"
#include "domain/mesh/mesh_i.h"
#include "problem/parameter_types.h"

namespace bart {

namespace domain {

/*! This struct provides the type of the triangulation to use. For 2D and 3D we
 * will use the built-in distributed triangulation, which handles distribution
 * of the triangulation across processors. For the 1D case, we need to make a
 * normal triangulation and use a separate function to distribute it. This is
 * a product of the way dealii is coded.
 */
template <int dim>
struct TriangulationType {
  using type = dealii::parallel::distributed::Triangulation<dim>;
};

template <>
struct TriangulationType<1> {
  using type = dealii::Triangulation<1>;
};

template struct TriangulationType<2>;
template struct TriangulationType<3>;

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
class Definition : public DefinitionI<dim> {
 public:
  typedef std::vector<typename dealii::DoFHandler<dim>::active_cell_iterator> CellRange;
  
  /*! \brief Constructor.
   * Takes ownership of injected dependencies (MeshI and FiniteElementI).
   */
  Definition(std::unique_ptr<domain::mesh::MeshI<dim>> mesh,
             std::shared_ptr<domain::finite_element::FiniteElementI<dim>> finite_element);
  ~Definition() = default;
  
  /*! Fills triangulation with mesh defined in MeshI object
   * Creates mesh shape, sets up boundary ids and material ids. Requires that 
   * the mesh has a material mapping setup.
   */ 
  Definition<dim>& SetUpMesh() override;

  /*! Set up the DOF handler, to access sparsity patterns, etc */
  Definition<dim>& SetUpDOF() override;

  /*! Get the parameters required to build a system matrix for this domain */
  void FillMatrixParameters(
      data::MatrixParameters &to_fill,
      problem::DiscretizationType discretization) const override;

  /*! Get a matrix suitible for a cell matrix.
   *
   * \return a dealii FullMatrix<double> of appropriate size.
   */
  dealii::FullMatrix<double> GetCellMatrix() const override {
    int cell_dofs = finite_element_->dofs_per_cell();
    dealii::FullMatrix<double> full_matrix(cell_dofs, cell_dofs);
    return full_matrix;
  }
  /*! Get a matrix suitible for a cell rhs.
   *
   * \return a dealii Vector<double> of appropriate size.
   */
  dealii::Vector<double> GetCellVector() const override {
    int cell_dofs = finite_element_->dofs_per_cell();
    dealii::Vector<double> vector(cell_dofs);
    return vector;
  }

  /*! Get a range of all cells to allow iterating over them */
  CellRange Cells() const override { return local_cells_; };

  /*! Get total degrees of freedom */
  int total_degrees_of_freedom() const override ;

  /*! Get locally owned degrees of freedom */
  dealii::IndexSet locally_owned_dofs() const {
    return locally_owned_dofs_;
  }

  /*! Get internal DOF object */
  const dealii::DoFHandler<dim>& dof_handler() const {
    return dof_handler_;
  }
  
 private:

  //! Internal owned mesh object.
  std::unique_ptr<domain::mesh::MeshI<dim>> mesh_;
  
  //! Internal owned finite element object
  std::shared_ptr<domain::finite_element::FiniteElementI<dim>> finite_element_;

  //! Internal distributed triangulation object
  typename TriangulationType<dim>::type triangulation_;

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

  /*! local cells */
  CellRange local_cells_;
};

} // namespace domain

} // namespace bart

#endif // BART_SRC_DOMAIN_DEFINITION_H_
