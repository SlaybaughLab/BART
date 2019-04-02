#include "../finite_element_gaussian.h"

#include <vector>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <gtest/gtest.h>

#include "problem/parameter_types.h"
#include "domain/tests/finite_element_test.h"
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

// BASE CLASS TESTS

class FiniteElementGaussianBaseMethods1D :
    public bart::domain::testing::FiniteElementBaseClassTest<1> {
 protected:
  using DiscretizationType = bart::problem::DiscretizationType;
  void SetUp() override {
    bart::domain::testing::FiniteElementBaseClassTest<1>::SetUp();
  }
};

TEST_F(FiniteElementGaussianBaseMethods1D, BaseTests) {
  bart::domain::FiniteElementGaussian<1> test_fe{DiscretizationType::kDiscontinuousFEM, 2};
  TestSetCell(&test_fe);
  TestSetCellAndFace(&test_fe);
}

class FiniteElementGaussianBaseMethods2D :
    public bart::domain::testing::FiniteElementBaseClassTest<2> {
 protected:
  using DiscretizationType = bart::problem::DiscretizationType;
  void SetUp() override {
    bart::domain::testing::FiniteElementBaseClassTest<2>::SetUp();
  }
};

TEST_F(FiniteElementGaussianBaseMethods2D, BaseTests) {
  bart::domain::FiniteElementGaussian<2> test_fe{DiscretizationType::kDiscontinuousFEM, 2};
  TestSetCell(&test_fe);
  TestSetCellAndFace(&test_fe);
}

TEST_F(FiniteElementGaussianTest, BasisValueTest) {
  bart::domain::FiniteElementGaussian<2> test_fe{DiscretizationType::kDiscontinuousFEM, 2};
  dealii::Triangulation<2> triangulation;

  dealii::GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(2);

  dealii::DoFHandler dof_handler(triangulation);
  dof_handler.distribute_dofs(*test_fe.finite_element());

  test_fe.values()->reinit(dof_handler.begin_active());

  int cell_dofs = test_fe.dofs_per_cell();
  int cell_quad_points = test_fe.n_cell_quad_pts();

  for (int i = 0; i < cell_dofs; ++i) {
    for (int q = 0; q < cell_quad_points; ++q) {
      EXPECT_DOUBLE_EQ(test_fe.values()->shape_value(i, q),
                       test_fe.ShapeValue(i, q));

    }
  }
}



} // namespace