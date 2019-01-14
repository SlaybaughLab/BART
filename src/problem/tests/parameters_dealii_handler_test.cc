#include "../parameters_dealii_handler.h"

#include <stdexcept>
#include <vector>

#include "gtest/gtest.h"

class ParametersDealiiHandlerTest : public ::testing::Test {
 protected:
  void SetUp() override;
  dealii::ParameterHandler test_parameter_handler;
  bart::problem::ParametersDealiiHandler test_parameters;

  // Key-words for input file
  const std::string kSpatialDimension_ = "problem dimension";
  const std::string kSpatialMax_ = "x, y, z max values of boundary locations";
};

void ParametersDealiiHandlerTest::SetUp() {
  test_parameters.SetUp(test_parameter_handler);
}

TEST_F(ParametersDealiiHandlerTest, BasicParametersDefault) {
  test_parameters.Parse(test_parameter_handler);
  
  ASSERT_EQ(test_parameters.SpatialDimension(), 2)
      << "Default spatial dimension";
}

TEST_F(ParametersDealiiHandlerTest, BasicParametersParse) {

  std::vector<double> spatial_max{10.0, 5.0, 8.0};
  
  test_parameter_handler.set(kSpatialDimension_, 3.0);
  test_parameter_handler.set(kSpatialMax_, "10.0, 5.0, 8.0");
  test_parameters.Parse(test_parameter_handler);

  ASSERT_EQ(test_parameters.SpatialDimension(), 3.0)
      << "Parsed spatial dimension";
  ASSERT_EQ(test_parameters.SpatialMax(), spatial_max)            
      << "Parsed spatial dimension";
}
