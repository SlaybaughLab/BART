#include "convergence/final_checker_or_n.h"

#include <memory>

#include <gtest/gtest.h>

#include "convergence/status.h"
#include "convergence/moments/tests/single_moment_checker_mock.h"
#include "data/moment_types.h"
#include "test_helpers/gmock_wrapper.h"

namespace {

using ::testing::_;
using ::testing::Return;
using namespace bart::convergence;

class ConvergenceFinalCheckerOrNSingleMomentTest : public ::testing::Test {
 protected:
  using FinalSingleMomentChecker =
      FinalCheckerOrN<bart::data::MomentVector, moments::SingleMomentCheckerI>;
  std::unique_ptr<moments::SingleMomentCheckerMock> checker_ptr;
  void SetUp() override;
};

void ConvergenceFinalCheckerOrNSingleMomentTest::SetUp() {
  checker_ptr = std::make_unique<moments::SingleMomentCheckerMock>();
  ON_CALL(*checker_ptr, CheckIfConverged(_,_)).
      WillByDefault(Return(true));
}

TEST_F(ConvergenceFinalCheckerOrNSingleMomentTest, Constructor) {
  FinalSingleMomentChecker test_checker(std::move(checker_ptr));

  EXPECT_EQ(checker_ptr, nullptr);
}

// Verify setters and getters work properly
TEST_F(ConvergenceFinalCheckerOrNSingleMomentTest, SettersAndGetters) {
  FinalSingleMomentChecker test_checker(std::move(checker_ptr));

  EXPECT_EQ(test_checker.max_iterations(), 100);
  test_checker.SetMaxIterations(50);
  EXPECT_EQ(test_checker.max_iterations(), 50);

  EXPECT_EQ(test_checker.iteration(), 0);
  test_checker.SetIteration(10);
  EXPECT_EQ(test_checker.iteration(), 10);
}

// -- CONVERGENCE TESTS --

// -- ERRORS --
// Verify setters and getters throw the correct errors
TEST_F(ConvergenceFinalCheckerOrNSingleMomentTest, SettersBadValues) {
  FinalSingleMomentChecker test_checker(std::move(checker_ptr));

  EXPECT_ANY_THROW(test_checker.SetIteration(-10));
  EXPECT_ANY_THROW(test_checker.SetMaxIterations(-10));
  EXPECT_ANY_THROW(test_checker.SetMaxIterations(0));

}


} // namespace




