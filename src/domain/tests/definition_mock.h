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

  MOCK_METHOD(DefinitionMock<dim>&, SetUpMesh, (), (override));
  MOCK_METHOD(DefinitionMock<dim>&, SetUpDOF, (), (override));
  MOCK_METHOD(void, FillMatrixParameters,
      (data::MatrixParameters &, problem::DiscretizationType),
      (override, const));
  MOCK_METHOD(dealii::FullMatrix<double>, GetCellMatrix, (), (override, const));
  MOCK_METHOD(dealii::Vector<double>, GetCellVector, (), (override, const));
  MOCK_METHOD(typename DefinitionI<dim>::CellRange, Cells, (), (override, const));
  MOCK_METHOD(problem::DiscretizationType, discretization_type, (), (override, const));
  MOCK_METHOD(int, total_degrees_of_freedom, (), (override, const));
  MOCK_METHOD(const dealii::DoFHandler<dim>&, dof_handler, (), (override, const));
  MOCK_METHOD(dealii::IndexSet, locally_owned_dofs, (), (override, const));

  };

} // namespace domain

} // namespace bart

#endif // BART_SRC_DOMAIN_TESTS_DEFINITION_MOCK_H_