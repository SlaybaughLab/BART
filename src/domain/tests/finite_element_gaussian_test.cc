#include "../finite_element_gaussian.h"

#include <vector>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <gtest/gtest.h>

#include "problem/parameter_types.h"

#include "test_helpers/gmock_wrapper.h"

namespace {

class FiniteElementGaussianTest : public ::testing::Test {
 protected:
  using DiscretizationType = bart::problem::DiscretizationType;
};

TEST_F(FiniteElementGaussianTest, ConstructorContinuous) {
  bart::domain::FiniteElementGaussian<2> test_fe{DiscretizationType::kContinuousFEM, 2};
  auto fe_q_ptr =
      dynamic_cast<dealii::FE_Q<2>*>(test_fe.finite_element());
  auto fe_value_ptr =
      dynamic_cast<dealii::FEValues<2>*>(test_fe.values());
  auto fe_face_value_ptr = dynamic_cast<dealii::FEFaceValues<2>*>(
      test_fe.face_values());
  auto fe_neighbor_face_value_ptr = dynamic_cast<dealii::FEFaceValues<2>*>(
      test_fe.neighbor_face_values());
  auto cell_quad_ptr =
      dynamic_cast<dealii::QGauss<2>*>(test_fe.cell_quadrature());
  auto face_quad_ptr =
      dynamic_cast<dealii::QGauss<1>*>(test_fe.face_quadrature());
  
  ASSERT_EQ(test_fe.polynomial_degree(), 2);
  ASSERT_FALSE(fe_q_ptr == nullptr);
  ASSERT_FALSE(fe_value_ptr == nullptr);
  ASSERT_FALSE(fe_face_value_ptr == nullptr);
  ASSERT_TRUE(fe_neighbor_face_value_ptr == nullptr);
  ASSERT_FALSE(cell_quad_ptr == nullptr);
  ASSERT_FALSE(face_quad_ptr == nullptr);
  ASSERT_EQ(test_fe.dofs_per_cell(), 9);
  ASSERT_EQ(test_fe.n_cell_quad_pts(), 9);
  ASSERT_EQ(test_fe.n_face_quad_pts(), 3);
}

TEST_F(FiniteElementGaussianTest, ConstructorDiscontinuous) {
  bart::domain::FiniteElementGaussian<2> test_fe{DiscretizationType::kDiscontinuousFEM, 2};
  dealii::FE_DGQ<2> *fe_q_ptr =
      dynamic_cast<dealii::FE_DGQ<2>*>(test_fe.finite_element());
  auto fe_neighbor_face_value_ptr = dynamic_cast<dealii::FEFaceValues<2>*>(
      test_fe.neighbor_face_values());
  ASSERT_FALSE(fe_q_ptr == nullptr);
  ASSERT_FALSE(fe_neighbor_face_value_ptr == nullptr);
}

TEST_F(FiniteElementGaussianTest, ConstructorNone) {
  ASSERT_ANY_THROW({
      bart::domain::FiniteElementGaussian<2> test_fe(DiscretizationType::kNone, 2);
    });
}

TEST_F(FiniteElementGaussianTest, SetCellTest) {
  bart::domain::FiniteElementGaussian<2> test_fe{DiscretizationType::kDiscontinuousFEM, 2};
  dealii::Triangulation<2> triangulation;

  dealii::GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(2);

  dealii::DoFHandler dof_handler(triangulation);
  dof_handler.distribute_dofs(*test_fe.finite_element());

  auto cell = dof_handler.begin_active();
  auto cell_id = cell->id();

  EXPECT_NO_THROW(test_fe.SetCell(cell));

  test_fe.values()->reinit(cell);

  EXPECT_FALSE(test_fe.SetCell(cell)); // Shouldn't change anything
  EXPECT_EQ(cell_id, test_fe.values()->get_cell()->id()); // Cell didn't change

  auto next_cell = cell;
  ++next_cell;
  auto next_cell_id = next_cell->id();

  EXPECT_TRUE(test_fe.SetCell(next_cell));
  // Check changed
  EXPECT_NE(cell_id, test_fe.values()->get_cell()->id());
  EXPECT_EQ(next_cell_id, test_fe.values()->get_cell()->id());
}

TEST_F(FiniteElementGaussianTest, SetCellAndFace) {
  bart::domain::FiniteElementGaussian<2> test_fe{DiscretizationType::kDiscontinuousFEM, 2};
  dealii::Triangulation<2> triangulation;

  dealii::GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(2);

  dealii::DoFHandler dof_handler(triangulation);
  dof_handler.distribute_dofs(*test_fe.finite_element());

  auto cell = dof_handler.begin_active();
  auto cell_id = cell->id();
  int face = 0;
  int face_index = cell->face_index(face);

  test_fe.face_values()->reinit(cell, face);

  EXPECT_FALSE(test_fe.SetFace(cell, face));
  EXPECT_EQ(cell_id, test_fe.face_values()->get_cell()->id());
  EXPECT_EQ(face_index, test_fe.face_values()->get_face_index());

  auto next_cell = cell;
  ++next_cell;
  auto next_cell_id = next_cell->id();
  int next_face = face + 1;
  int next_face_index = next_cell->face_index(next_face);

  EXPECT_TRUE(test_fe.SetFace(next_cell, next_face));
  EXPECT_EQ(next_cell_id, test_fe.face_values()->get_cell()->id());
  EXPECT_NE(face_index, test_fe.face_values()->get_face_index());
  EXPECT_EQ(next_face_index, test_fe.face_values()->get_face_index());
}



} // namespace