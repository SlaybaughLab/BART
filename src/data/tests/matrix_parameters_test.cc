#include "data/matrix_parameters.h"

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>

#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

class MatrixParametersTest : public ::testing::Test {

};

TEST_F(MatrixParametersTest, BuildMatrixTest) {
  dealii::Triangulation<2> triangulation;
  dealii::FE_Q<2> finite_element(1);
  dealii::DoFHandler<2> dof_handler(triangulation);

  dealii::GridGenerator::hyper_cube(triangulation, 0, 1);
  dof_handler.distribute_dofs(finite_element);

  data::MatrixParameters test_parameters;
  test_parameters.sparsity_pattern.reinit(dof_handler.n_dofs(), dof_handler.n_dofs());
  dealii::DoFTools::make_sparsity_pattern(dof_handler, test_parameters.sparsity_pattern);
  test_parameters.rows = dof_handler.locally_owned_dofs();
  test_parameters.columns = dof_handler.locally_owned_dofs();

  auto matrix_pointer = data::BuildMatrix(test_parameters);

  EXPECT_EQ(matrix_pointer->m(), dof_handler.n_dofs());
  EXPECT_EQ(matrix_pointer->n(), dof_handler.n_dofs());
}

} // namespace