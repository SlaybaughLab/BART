#include "convergence/final_checker_or_n.h"

#include <memory>
#include <optional>

#include <gtest/gtest.h>

#include "convergence/status.h"
#include "convergence/reporter/tests/mpi_mock.h"
#include "convergence/moments/tests/multi_moment_checker_mock.h"
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

class ConvergenceFinalCheckerOrNMultiMomentTest :
    public bart::convergence::testing::ConvergenceFinalTest<bart::system::moments::MomentsMap> {
 protected:
  using Reporter = reporter::MpiMock;
  using FinalMultiMomentChecker =
      FinalCheckerOrN<bart::system::moments::MomentsMap , moments::MultiMomentCheckerI>;
  
  std::unique_ptr<NiceMock<moments::MultiMomentCheckerMock>> checker_ptr;
  std::shared_ptr<Reporter> reporter_ptr;

  bart::system::moments::MomentsMap moment_map_one, moment_map_two;
  void SetUp() override;
};

void ConvergenceFinalCheckerOrNMultiMomentTest::SetUp() {
  checker_ptr = std::make_unique<NiceMock<moments::MultiMomentCheckerMock>>();
  reporter_ptr = std::make_shared<Reporter>();
  ON_CALL(*checker_ptr, CheckIfConverged(_,_))
      .WillByDefault(Return(true));
  ON_CALL(*checker_ptr, delta())
      .WillByDefault(Return(std::nullopt));
  ON_CALL(*checker_ptr, failed_index())
      .WillByDefault(Return(std::nullopt));
}

TEST_F(ConvergenceFinalCheckerOrNMultiMomentTest, BaseClassTests) {
  FinalMultiMomentChecker test_checker(std::move(checker_ptr));
  TestBaseMethods(&test_checker);
}

TEST_F(ConvergenceFinalCheckerOrNMultiMomentTest, Constructor) {
  FinalMultiMomentChecker test_checker(std::move(checker_ptr));

  EXPECT_EQ(checker_ptr, nullptr);
  EXPECT_NE(nullptr,
            dynamic_cast<NiceMock<moments::MultiMomentCheckerMock>*>(
                test_checker.checker_ptr()));
}

// -- CONVERGENCE TESTS --

TEST_F(ConvergenceFinalCheckerOrNMultiMomentTest, GoodConvergence) {
  FinalMultiMomentChecker test_checker(std::move(checker_ptr));
  Status good_convergence = {1, 100, true, std::nullopt, std::nullopt};

  auto result = test_checker.CheckFinalConvergence(moment_map_one, moment_map_two);
  EXPECT_TRUE(CompareStatus(result, good_convergence));
  EXPECT_TRUE(CompareStatus(test_checker.convergence_status(),
                            good_convergence));
  EXPECT_TRUE(test_checker.convergence_is_complete());
}

TEST_F(ConvergenceFinalCheckerOrNMultiMomentTest, GoodConvergenceAfterBad) {
  Expectation bad_convergence = EXPECT_CALL(*checker_ptr, CheckIfConverged(_,_))
      .Times(5)
      .WillRepeatedly(Return(false));
  EXPECT_CALL(*checker_ptr, CheckIfConverged(_,_))
      .After(bad_convergence)
      .WillOnce(Return(true));
  EXPECT_CALL(*reporter_ptr, Report(A<const Status&>()))
      .Times(6);

  FinalMultiMomentChecker test_checker(std::move(checker_ptr), reporter_ptr);
  Status result, good_convergence = {6, 100, true, std::nullopt, std::nullopt};

  for (int i = 0; i < 6; ++i)
    result = test_checker.CheckFinalConvergence(moment_map_one, moment_map_two);

  EXPECT_TRUE(CompareStatus(result, good_convergence));
  EXPECT_TRUE(test_checker.convergence_is_complete());
}

TEST_F(ConvergenceFinalCheckerOrNMultiMomentTest, BadConvergenceAfterGood) {

  auto delta = std::make_optional<double>(0.123);
  auto failed_index = std::make_optional<int>(5);

  Expectation good_convergence = EXPECT_CALL(*checker_ptr, CheckIfConverged(_,_))
      .Times(5)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*checker_ptr, delta())
      .Times(5)
      .WillRepeatedly(::testing::DoDefault());
  EXPECT_CALL(*checker_ptr, failed_index())
      .Times(5)
      .WillRepeatedly(::testing::DoDefault());
  Expectation bad_convergence = EXPECT_CALL(*checker_ptr, CheckIfConverged(_,_))
      .After(good_convergence)
      .WillOnce(Return(false));
  EXPECT_CALL(*checker_ptr, delta())
      .After(bad_convergence)
      .WillOnce(Return(delta));
  EXPECT_CALL(*checker_ptr, failed_index())
      .After(bad_convergence)
      .WillOnce(Return(failed_index));
  EXPECT_CALL(*reporter_ptr, Report(A<const Status&>()))
      .Times(6);

  FinalMultiMomentChecker test_checker(std::move(checker_ptr), reporter_ptr);
  Status result, expected = {6, 100, false, failed_index, delta};

  for (int i = 0; i < 6; ++i)
    result = test_checker.CheckFinalConvergence(moment_map_one, moment_map_two);

  EXPECT_TRUE(CompareStatus(result, expected));
  EXPECT_FALSE(test_checker.convergence_is_complete());
}
//
TEST_F(ConvergenceFinalCheckerOrNMultiMomentTest, MaxIterationsReached) {

  auto delta = std::make_optional<double>(0.123);
  auto failed_index = std::make_optional<int>(5);

  ON_CALL(*checker_ptr, CheckIfConverged(_,_))
      .WillByDefault(Return(false));
  ON_CALL(*checker_ptr, delta())
      .WillByDefault(Return(delta));
  ON_CALL(*checker_ptr, failed_index())
      .WillByDefault(Return(failed_index));

  FinalMultiMomentChecker test_checker(std::move(checker_ptr));
  Status result, expected = {10, 10, true, failed_index, delta};

  test_checker.SetMaxIterations(10).SetIteration(9);

  result = test_checker.CheckFinalConvergence(moment_map_one, moment_map_two);

  EXPECT_TRUE(CompareStatus(result, expected));
  EXPECT_TRUE(test_checker.convergence_is_complete());
}

} // namespace




