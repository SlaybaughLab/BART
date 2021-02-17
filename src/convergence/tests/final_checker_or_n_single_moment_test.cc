#include "convergence/final_checker_or_n.h"

#include <memory>
#include <optional>

#include <gtest/gtest.h>

#include "convergence/status.hpp"
#include "convergence/moments/tests/single_moment_checker_mock.h"
#include "convergence/tests/final_test.h"
#include "system/moments/spherical_harmonic_types.h"
#include "test_helpers/gmock_wrapper.h"

namespace {

using ::testing::_;
using ::testing::Expectation;
using ::testing::NiceMock;
using ::testing::Return;
using ::testing::A;
using namespace bart::convergence;
using bart::convergence::testing::CompareStatus;
using bart::system::moments::MomentVector;


class ConvergenceFinalCheckerOrNSingleMomentTest :
    public bart::convergence::testing::ConvergenceFinalTest<MomentVector> {
 protected:
  using FinalSingleMomentChecker =
      FinalCheckerOrN<MomentVector, moments::SingleMomentCheckerI>;
  std::unique_ptr<NiceMock<moments::SingleMomentCheckerMock>> checker_ptr;
  bart::system::moments::MomentVector moment_one, moment_two;
  void SetUp() override;
};

void ConvergenceFinalCheckerOrNSingleMomentTest::SetUp() {
  checker_ptr = std::make_unique<NiceMock<moments::SingleMomentCheckerMock>>();
  ON_CALL(*checker_ptr, IsConverged(_,_))
      .WillByDefault(Return(true));
  ON_CALL(*checker_ptr, delta())
      .WillByDefault(Return(std::nullopt));
}

// -- BASE CLASS METHODS

TEST_F(ConvergenceFinalCheckerOrNSingleMomentTest, BaseClass) {
  FinalSingleMomentChecker test_checker(std::move(checker_ptr));
  TestBaseMethods(&test_checker);
}

TEST_F(ConvergenceFinalCheckerOrNSingleMomentTest, Constructor) {
  FinalSingleMomentChecker test_checker(std::move(checker_ptr));

  EXPECT_EQ(checker_ptr, nullptr);
  EXPECT_NE(nullptr,
            dynamic_cast<NiceMock<moments::SingleMomentCheckerMock>*>(
                test_checker.checker_ptr()));
}

// -- CONVERGENCE TESTS --

TEST_F(ConvergenceFinalCheckerOrNSingleMomentTest, GoodConvergence) {
  FinalSingleMomentChecker test_checker(std::move(checker_ptr));
  Status good_convergence = {1, 100, true, std::nullopt, std::nullopt};
  auto result = test_checker.ConvergenceStatus(moment_one, moment_two);
  EXPECT_TRUE(CompareStatus(result, good_convergence));
}

TEST_F(ConvergenceFinalCheckerOrNSingleMomentTest, GoodConvergenceAfterBad) {
  Expectation bad_convergence = EXPECT_CALL(*checker_ptr, IsConverged(_,_))
      .Times(5)
      .WillRepeatedly(Return(false));
  EXPECT_CALL(*checker_ptr, IsConverged(_,_))
      .After(bad_convergence)
      .WillOnce(Return(true));

  FinalSingleMomentChecker test_checker(std::move(checker_ptr));
  Status result, good_convergence = {6, 100, true, std::nullopt, std::nullopt};

  for (int i = 0; i < 6; ++i)
    result = test_checker.ConvergenceStatus(moment_one, moment_two);

  EXPECT_TRUE(CompareStatus(result, good_convergence));
}

TEST_F(ConvergenceFinalCheckerOrNSingleMomentTest, BadConvergenceAfterGood) {

  auto delta = std::make_optional<double>(0.123);

  Expectation good_convergence = EXPECT_CALL(*checker_ptr, IsConverged(_,_))
      .Times(5)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*checker_ptr, delta())
      .Times(5)
      .WillRepeatedly(::testing::DoDefault());
  Expectation bad_convergence = EXPECT_CALL(*checker_ptr, IsConverged(_,_))
      .After(good_convergence)
      .WillOnce(Return(false));
  EXPECT_CALL(*checker_ptr, delta())
      .After(bad_convergence)
      .WillOnce(Return(delta));

  FinalSingleMomentChecker test_checker(std::move(checker_ptr));
  Status result, expected = {6, 100, false, std::nullopt, delta};

  for (int i = 0; i < 6; ++i)
    result = test_checker.ConvergenceStatus(moment_one, moment_two);

  EXPECT_TRUE(CompareStatus(result, expected));
}

TEST_F(ConvergenceFinalCheckerOrNSingleMomentTest, MaxIterationsReached) {

  auto delta = std::make_optional<double>(0.123);

  ON_CALL(*checker_ptr, IsConverged(_,_))
      .WillByDefault(Return(false));
  ON_CALL(*checker_ptr, delta())
      .WillByDefault(Return(delta));

  FinalSingleMomentChecker test_checker(std::move(checker_ptr));
  Status result, expected = {10, 10, true, std::nullopt, delta};

  test_checker.SetMaxIterations(10).SetIteration(9);

  result = test_checker.ConvergenceStatus(moment_one, moment_two);

  EXPECT_TRUE(CompareStatus(result, expected));

  Status reset;
  reset.max_iterations = expected.max_iterations;
  test_checker.Reset();
  EXPECT_TRUE(CompareStatus(test_checker.convergence_status(), reset));
}


} // namespace




