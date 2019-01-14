#include "../parameters_dealii_handler.h"

#include "gtest/gtest.h"

class ParametersDealiiHandlerTest : public ::testing::Test {
 protected:
  void SetUp() override;
  dealii::ParameterHandler test_parameter_handler;
  bart::problem::ParametersDealiiHandler test_parameters;
};

void ParametersDealiiHandlerTest::SetUp() {
  test_parameters.SetUp(test_parameter_handler);
}

TEST_F(ParametersDealiiHandlerTest, BasicParametersParsing) {
  test_parameters.Parse(test_parameter_handler);
  
  ASSERT_EQ(test_parameters.SpatialDimension(), 2);
}
