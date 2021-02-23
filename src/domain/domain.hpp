#ifndef BART_SRC_DOMAIN_DOMAIN_HPP_
#define BART_SRC_DOMAIN_DOMAIN_HPP_

#include <memory>
#include <vector>

#include <deal.II/base/index_set.h>
#include <deal.II/base/iterator_range.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include "domain/domain_i.hpp"
#include "domain/finite_element/finite_element_i.hpp"
#include "domain/mesh/mesh_i.h"
#include "problem/parameter_types.hpp"
#include "system/system_types.h"

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

/*! \brief Default implementation of a domain definition.
 *
 * See the base class for general information about domain definitions, an
 * example of use is shown below:
 *
 * \code{.cpp}
 *
 * domain::definition<3> domain_definition(std::move(mesh_ptr),
 *                                         finite_element_ptr,
 *                                         problem::DiscretizationType::kContinuousFEM);
 * domain_definition.SetUpMesh().SetUpDOF();
 *
 * auto cell_matrix = domain_definition.GetCellMatrix();
 *
 * \endcode
 *
 * \author Joshua Rehak
 * \date 2019/02
 */
template <int dim>
class Definition : public DefinitionI<dim> {
 public:
  typedef std::vector<typename dealii::DoFHandler<dim>::active_cell_iterator> CellRange;
  
  /*! \brief Constructor.
   * Takes ownership of injected dependencies (MeshI and FiniteElementI) and
   * sets the type of discretization (default: continuous FEM).
   */
  Definition(std::unique_ptr<domain::mesh::MeshI<dim>> mesh,
             std::shared_ptr<domain::finite_element::FiniteElementI<dim>> finite_element,
             problem::DiscretizationType discretization = problem::DiscretizationType::kContinuousFEM);
  ~Definition() = default;

  Definition<dim>& SetUpDOF() override;
  Definition<dim>& SetUpMesh() override;
  Definition<dim>& SetUpMesh(const int global_refinements) override;

  dealii::FullMatrix<double> GetCellMatrix() const override {
    int cell_dofs = finite_element_->dofs_per_cell();
    dealii::FullMatrix<double> full_matrix(cell_dofs, cell_dofs);
    return full_matrix;
  }

  dealii::Vector<double> GetCellVector() const override {
    int cell_dofs = finite_element_->dofs_per_cell();
    dealii::Vector<double> vector(cell_dofs);
    return vector;
  }

  std::shared_ptr<system::MPISparseMatrix> MakeSystemMatrix() const override;

  std::shared_ptr<system::MPIVector> MakeSystemVector() const override;

  CellRange Cells() const override { return local_cells_; };

  problem::DiscretizationType discretization_type() const override {
    return discretization_type_; }

  int total_degrees_of_freedom() const override ;

  dealii::IndexSet locally_owned_dofs() const override {
    return locally_owned_dofs_; }


  const dealii::DoFHandler<dim>& dof_handler() const override {
    return dof_handler_; }
  
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
  dealii::AffineConstraints<double> constraint_matrix_;

  /*! Dynamic sparsity pattern for MPI matrices */
  dealii::DynamicSparsityPattern dynamic_sparsity_pattern_;

  /*! local cells */
  CellRange local_cells_;

  /*! Discretization type */
  const problem::DiscretizationType discretization_type_;
};

} // namespace domain

} // namespace bart

#endif // BART_SRC_DOMAIN_DOMAIN_HPP_
