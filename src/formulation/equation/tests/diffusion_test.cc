#include "formulation/types.h"
#include "formulation/equation/diffusion.h"

#include <gtest/gtest.h>

#include "test_helpers/gmock_wrapper.h"

class EquationDiffusionTest : public ::testing::Test {
 protected:
  using DiscretizationType = bart::formulation::DiscretizationType;
  using ScalarEquations = bart::formulation::ScalarEquations;
  using EquationType = bart::formulation::EquationType;

  bart::formulation::equation::Diffusion<2> diffusion_cfem{
    DiscretizationType::kContinuousFEM};
  bart::formulation::equation::Diffusion<2> diffusion_dfem{
    DiscretizationType::kDiscontinuousFEM};
};

TEST_F(EquationDiffusionTest, ConstructorTest) {
  EXPECT_EQ(diffusion_cfem.discretization_type(),
      DiscretizationType::kContinuousFEM);
  EXPECT_EQ(diffusion_cfem.equation_type(), EquationType::kScalar);

  EXPECT_EQ(diffusion_dfem.discretization_type(),
            DiscretizationType::kDiscontinuousFEM);
  EXPECT_EQ(diffusion_dfem.equation_type(), EquationType::kScalar);
}