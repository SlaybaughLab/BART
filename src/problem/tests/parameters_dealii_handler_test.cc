#include "../parameters_dealii_handler.h"

#include "gtest/gtest.h"

class ParametersDealiiHandlerTest : public ::testing::Test {
 protected:
  void SetUp() override;
  dealii::ParameterHandler test_prm;
};

void ParametersDealiiHandlerTest::SetUp() {
  test_prm.declare_entry ("problem dimension", "2",
                          dealii::Patterns::Integer(), "");
}

TEST_F(ParametersDealiiHandlerTest, BasicParametersParsing) {
  bart::problem::ParametersDealiiHandler test_parameters{test_prm};
  ASSERT_EQ(test_parameters.SpatialDimension(), 2);
}
