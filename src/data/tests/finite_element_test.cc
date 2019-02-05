#include "../finite_element.h"

#include <gtest/gtest.h>

#include "../../test_helpers/gmock_wrapper.h"

class FiniteElementTest : public ::testing::Test {
 protected:
  using DiscretizationType = bart::problem::DiscretizationType;
};

TEST_F(FiniteElementTest, ConstructorPolynomialDegree) {
  bart::data::FiniteElement<1> test_fe{DiscretizationType::kContinuousFEM, 2};
  ASSERT_EQ(test_fe.polynomial_degree(), 2);
}
