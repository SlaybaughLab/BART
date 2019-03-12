#include "problem/parameter_types.h"
#include "formulation/types.h"
#include "formulation/equation/diffusion.h"

#include <gtest/gtest.h>

#include "test_helpers/gmock_wrapper.h"

class EquationDiffusionTest : public ::testing::Test {
 protected:
  using DiscretizationType = bart::formulation::DiscretizationType;
  using ScalarEquations = bart::formulation::ScalarEquations;
  using EquationType = bart::formulation::EquationType;
  using ProblemType = bart::problem::ProblemType;

};

TEST_F(EquationDiffusionTest, ConstructorTest) {
  auto k_effective_ptr = std::make_shared<double>(1.0);

  bart::formulation::equation::Diffusion<2> diffusion_cfem{
      DiscretizationType::kContinuousFEM,
      ProblemType::kFixedSource,
      k_effective_ptr};

  bart::formulation::equation::Diffusion<2> diffusion_dfem{
      DiscretizationType::kDiscontinuousFEM,
      ProblemType::kFixedSource,
      k_effective_ptr};

  EXPECT_EQ(diffusion_cfem.discretization_type(),
      DiscretizationType::kContinuousFEM);
  EXPECT_EQ(diffusion_cfem.equation_type(), EquationType::kScalar);

  EXPECT_EQ(diffusion_dfem.discretization_type(),
            DiscretizationType::kDiscontinuousFEM);
  EXPECT_EQ(diffusion_dfem.equation_type(), EquationType::kScalar);
}