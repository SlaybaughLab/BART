
#include "convergence/checker_factory.h"

#include "system/moments/spherical_harmonic_types.h"
#include "convergence/moments/single_moment_checker_l1_norm.h"
#include "convergence/final_checker_or_n.h"
#include "convergence/parameters/single_parameter_checker.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class ConvergenceCheckerFactoryIntegrationTest : public ::testing::Test {
 public:
};

TEST_F(ConvergenceCheckerFactoryIntegrationTest, ParameterChecker) {
  const double max_delta = 1e-10;
  const int max_iterations = 20;

  auto checker_ptr = convergence::CheckerFactory::MakeParameterChecker(
          max_delta, max_iterations);

  using ExpectedType = convergence::FinalCheckerOrN<
      double,
      convergence::parameters::SingleParameterChecker>;

  ASSERT_NE(nullptr, checker_ptr);
  auto final_checker_ptr =
      dynamic_cast<ExpectedType*>(checker_ptr.get());
  ASSERT_NE(nullptr, final_checker_ptr);
  EXPECT_EQ(final_checker_ptr->max_iterations(), max_iterations);
  EXPECT_EQ(final_checker_ptr->checker_ptr()->max_delta(), max_delta);
}

TEST_F(ConvergenceCheckerFactoryIntegrationTest, SingleMomentChecker) {
  const double max_delta = 1e-10;
  const int max_iterations = 20;

  auto checker_ptr = convergence::CheckerFactory::MakeSingleMomentChecker(
          max_delta, max_iterations);

  using ExpectedType = convergence::FinalCheckerOrN<
      system::moments::MomentVector,
      convergence::moments::SingleMomentCheckerI>;

  ASSERT_NE(nullptr, checker_ptr);
  auto final_checker_ptr =
      dynamic_cast<ExpectedType*>(checker_ptr.get());
  ASSERT_NE(nullptr, final_checker_ptr);
  EXPECT_EQ(final_checker_ptr->max_iterations(), max_iterations);

  using ExpectedCheckerType = convergence::moments::SingleMomentCheckerL1Norm;
  auto single_checker_ptr =
      dynamic_cast<ExpectedCheckerType*>(final_checker_ptr->checker_ptr());

  ASSERT_NE(nullptr, single_checker_ptr);
  EXPECT_EQ(single_checker_ptr->max_delta(), max_delta);
}


} // namespace