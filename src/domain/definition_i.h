#ifndef BART_SRC_DOMAIN_DEFINITION_I_H_
#define BART_SRC_DOMAIN_DEFINITION_I_H_

#include <vector>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/lac/full_matrix.h>

#include "data/matrix_parameters.h"
#include "problem/parameter_types.h"

namespace bart {

namespace domain {

template <int dim>
class DefinitionI {
 public:

  using Cell = typename dealii::DoFHandler<dim>::active_cell_iterator;
  using CellRange = std::vector<Cell>;

  virtual ~DefinitionI() = default;

  virtual DefinitionI<dim>& SetUpMesh() = 0;

  virtual DefinitionI<dim>& SetUpDOF() = 0;

  virtual void FillMatrixParameters(
      data::MatrixParameters &to_fill,
      problem::DiscretizationType discretization) const = 0;

  virtual dealii::FullMatrix<double> GetCellMatrix() const = 0;

  virtual CellRange Cells() const = 0;

  virtual int total_degrees_of_freedom() const = 0;
};

} // namespace domain

} // namespace bart

#endif // BART_SRC_DOMAIN_DEFINITION_I_H_