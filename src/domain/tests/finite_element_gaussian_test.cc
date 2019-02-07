#include "../finite_element_gaussian.h"

#include <vector>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <gtest/gtest.h>

#include "../../problem/parameter_types.h"

#include "../../test_helpers/gmock_wrapper.h"

class FiniteElementGaussianTest : public ::testing::Test {
 protected:
  using DiscretizationType = bart::problem::DiscretizationType;
};

TEST_F(FiniteElementGaussianTest, ConstructorContinuous) {
  bart::data::FiniteElementGaussian<2> test_fe{DiscretizationType::kContinuousFEM, 2};
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
  bart::data::FiniteElementGaussian<2> test_fe{DiscretizationType::kDiscontinuousFEM, 2};
  dealii::FE_DGQ<2> *fe_q_ptr =
      dynamic_cast<dealii::FE_DGQ<2>*>(test_fe.finite_element());
  auto fe_neighbor_face_value_ptr = dynamic_cast<dealii::FEFaceValues<2>*>(
      test_fe.neighbor_face_values());
  ASSERT_FALSE(fe_q_ptr == nullptr);
  ASSERT_FALSE(fe_neighbor_face_value_ptr == nullptr);
}

TEST_F(FiniteElementGaussianTest, ConstructorNone) {
  ASSERT_ANY_THROW({
      bart::data::FiniteElementGaussian<2> test_fe(DiscretizationType::kNone, 2);
    });
}
