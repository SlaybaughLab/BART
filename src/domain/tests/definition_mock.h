#ifndef BART_SRC_DOMAIN_TESTS_DEFINITION_MOCK_H_
#define BART_SRC_DOMAIN_TESTS_DEFINITION_MOCK_H_

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include "data/matrix_parameters.h"
#include "domain/definition_i.h"
#include "problem/parameter_types.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace domain {

template <int dim>
class DefinitionMock : public DefinitionI<dim> {
 public:
  using typename DefinitionI<dim>::CellRange;

  MOCK_METHOD0_T(SetUpMesh, DefinitionMock<dim>&());
  MOCK_METHOD0_T(SetUpDOF, DefinitionMock<dim>&());
  MOCK_CONST_METHOD2_T(FillMatrixParameters, void(data::MatrixParameters &to_fill,
      problem::DiscretizationType discretization));
  MOCK_CONST_METHOD0_T(GetCellMatrix, dealii::FullMatrix<double>());
  MOCK_CONST_METHOD0_T(GetCellVector, dealii::Vector<double>());
  MOCK_CONST_METHOD0_T(Cells, typename DefinitionI<dim>::CellRange());
  MOCK_CONST_METHOD0_T(total_degrees_of_freedom, int());
  MOCK_CONST_METHOD0_T(dof_handler, const dealii::DoFHandler<dim>&());

  };

} // namespace domain

} // namespace bart

#endif // BART_SRC_DOMAIN_TESTS_DEFINITION_MOCK_H_