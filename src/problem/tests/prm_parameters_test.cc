#include "../prm_parameters.h"

#include "gtest/gtest.h"

class PrmParametersTest : public ::testing::Test {
 protected:
  void SetUp() override;
  dealii::ParameterHandler test_prm;
};

void PrmParametersTest::SetUp() {
  test_prm.declare_entry ("problem dimension", "2",
                          dealii::Patterns::Integer(), "");
}

TEST_F(PrmParametersTest, BasicParametersParsing) {
  bart::problem::PrmParameters test_parameters{test_prm};
  ASSERT_EQ(test_parameters.Dimension(), 2);
}
