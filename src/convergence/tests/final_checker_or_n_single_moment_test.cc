#include "convergence/final_checker_or_n.h"

#include <memory>

#include <gtest/gtest.h>

#include "convergence/status.h"
#include "convergence/moments/tests/single_moment_checker_mock.h"
#include "convergence/moments/single_moment_checker_i.h"
#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart::convergence;

class ConvergenceFinalCheckerOrNSingleMomentTest : public ::testing::Test {
 protected:
  std::unique_ptr<moments::SingleMomentCheckerI> checker_ptr;
  void SetUp() override;
};

void ConvergenceFinalCheckerOrNSingleMomentTest::SetUp() {
  checker_ptr = std::make_unique<moments::SingleMomentCheckerMock>();
}

TEST_F(ConvergenceFinalCheckerOrNSingleMomentTest, Constructor) {
  FinalCheckerOrN<moments::SingleMomentCheckerI>
      test_checker(std::move(checker_ptr));

  EXPECT_EQ(checker_ptr, nullptr);
}

TEST_F(ConvergenceFinalCheckerOrNSingleMomentTest, SettersAndGetters) {
  FinalCheckerOrN<moments::SingleMomentCheckerI>
      test_checker(std::move(checker_ptr));

  EXPECT_EQ(test_checker.max_iterations(), 100);
  test_checker.SetMaxIterations(50);
  EXPECT_EQ(test_checker.max_iterations(), 50);

  EXPECT_EQ(test_checker.iteration(), 0);
  test_checker.SetIteration(10);
  EXPECT_EQ(test_checker.iteration(), 10);
}


} // namespace




