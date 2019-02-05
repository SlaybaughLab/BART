#include "../finite_element.h"

#include <vector>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <gtest/gtest.h>

#include "../../test_helpers/gmock_wrapper.h"

class FiniteElementTest : public ::testing::Test {
 protected:
  using DiscretizationType = bart::problem::DiscretizationType;
};

TEST_F(FiniteElementTest, ConstructorContinuous) {
  bart::data::FiniteElement<2> test_fe{DiscretizationType::kContinuousFEM, 2};
  dealii::FE_Q<2> *fe_q_ptr =
      dynamic_cast<dealii::FE_Q<2>*>(test_fe.finite_element());

  auto cell_quad_ptr =
      dynamic_cast<dealii::QGauss<2>*>(test_fe.cell_quadrature());
  auto face_quad_ptr =
      dynamic_cast<dealii::QGauss<1>*>(test_fe.face_quadrature());
  
  ASSERT_EQ(test_fe.polynomial_degree(), 2);
  ASSERT_FALSE(fe_q_ptr == nullptr);
  ASSERT_FALSE(cell_quad_ptr == nullptr);
  ASSERT_FALSE(face_quad_ptr == nullptr);  
}

TEST_F(FiniteElementTest, ConstructorDiscontinuous) {
  bart::data::FiniteElement<2> test_fe{DiscretizationType::kDiscontinuousFEM, 2};
  dealii::FE_DGQ<2> *fe_q_ptr =
      dynamic_cast<dealii::FE_DGQ<2>*>(test_fe.finite_element());

  auto cell_quad_ptr =
      dynamic_cast<dealii::QGauss<2>*>(test_fe.cell_quadrature());
  auto face_quad_ptr =
      dynamic_cast<dealii::QGauss<1>*>(test_fe.face_quadrature());
  
  ASSERT_EQ(test_fe.polynomial_degree(), 2);
  ASSERT_FALSE(fe_q_ptr == nullptr);
  ASSERT_FALSE(cell_quad_ptr == nullptr);
  ASSERT_FALSE(face_quad_ptr == nullptr);  
}

TEST_F(FiniteElementTest, ConstructorNone) {
  ASSERT_ANY_THROW({
      bart::data::FiniteElement<2> test_fe(DiscretizationType::kNone, 2);
    });
}
