#include "convergence/moments/single_moment_checker_l1_norm.h"


#include <gtest/gtest.h>

#include "convergence/tests/single_checker_test.h"
#include "data/moment_types.h"
#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"

namespace {

class SingleMomentCheckerL1NormTest :
    public bart::convergence::testing::SingleCheckerTest<bart::data::MomentVector> {
 protected:
  bart::convergence::moments::SingleMomentCheckerL1Norm checker{1e-6};
  bart::data::MomentVector moment_one;
  bart::data::MomentVector moment_two;
  void SetUp() override;
};

void SingleMomentCheckerL1NormTest::SetUp() {
  moment_one.reinit(5);
  moment_two.reinit(5);
  auto random_vector = btest::RandomVector(5, 0, 2);
  for (unsigned int i = 0; i < moment_one.size(); ++i) {
    moment_one(i) = random_vector[i];
    moment_two(i) = random_vector[i];
  }
}

TEST_F(SingleMomentCheckerL1NormTest, BaseMethods) {
  TestBaseMethods(&checker);
}

TEST_F(SingleMomentCheckerL1NormTest, SameVector) {
  EXPECT_TRUE(checker.CheckIfConverged(moment_one, moment_one));
  EXPECT_TRUE(checker.is_converged());
}

TEST_F(SingleMomentCheckerL1NormTest, OneThresholdAway) {
  double to_add = moment_one.l1_norm() * 0.99 * checker.max_delta();
  moment_two(2) += to_add;

  EXPECT_TRUE(checker.CheckIfConverged(moment_one, moment_two));
  EXPECT_TRUE(checker.CheckIfConverged(moment_two, moment_one));
  EXPECT_TRUE(checker.is_converged());
  EXPECT_NEAR(0.99 * checker.max_delta(),
              checker.delta().value(),
              1e-6);
}

TEST_F(SingleMomentCheckerL1NormTest, TwoThresholdAway) {
  double to_add = moment_one.l1_norm() * 2 * checker.max_delta();
  moment_two(2) += to_add;

  EXPECT_FALSE(checker.CheckIfConverged(moment_one, moment_two));
  EXPECT_FALSE(checker.CheckIfConverged(moment_two, moment_one));
  EXPECT_FALSE(checker.is_converged());
  EXPECT_NEAR(2 * checker.max_delta(),
              checker.delta().value(),
              1e-6);
}

TEST_F(SingleMomentCheckerL1NormTest, SetMaxDelta) {
  double to_set = 1e-5;
  checker.SetMaxDelta(to_set);
  EXPECT_EQ(checker.max_delta(), to_set);

  double to_add = moment_one.l1_norm() * 0.99 * to_set;
  moment_two(2) += to_add;

  EXPECT_TRUE(checker.CheckIfConverged(moment_one, moment_two));
  EXPECT_TRUE(checker.CheckIfConverged(moment_two, moment_one));

}

} // namespace