#include "convergence/final_flux_or_n.h"

#include <memory>
#include <optional>

#include "convergence/status.h"
#include "convergence/flux/tests/multi_checker_mock.h"
#include "data/system_fluxes.h"
#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"

using ::testing::_;
using ::testing::NiceMock;

class ConvergenceFinalFluxOrNTest : public ::testing::Test {
 protected:
  std::unique_ptr<bart::convergence::flux::MultiCheckerMock>
      mock_multi_checker_ptr;
  std::unique_ptr<NiceMock<bart::convergence::flux::MultiCheckerMock>>
      nice_mock_multi_checker_ptr;
  std::shared_ptr<bart::data::SystemFluxes> system_fluxes;
  void SetUp();
};

void ConvergenceFinalFluxOrNTest::SetUp() {
  mock_multi_checker_ptr =
      std::make_unique<bart::convergence::flux::MultiCheckerMock>();
  nice_mock_multi_checker_ptr =
      std::make_unique<NiceMock<bart::convergence::flux::MultiCheckerMock>>();
  system_fluxes = std::make_shared<bart::data::SystemFluxes>();
}    

TEST_F(ConvergenceFinalFluxOrNTest, DefaultStatus) {
  bart::convergence::FinalFluxOrN flux_tester(std::move(mock_multi_checker_ptr),
                                              system_fluxes);
  bart::convergence::Status status;
  status = flux_tester.convergence_status();
  ASSERT_EQ(status.iteration_number, 0);
  ASSERT_EQ(status.max_iterations, 1);
  ASSERT_FALSE(status.is_complete);
  ASSERT_FALSE(status.failed_index.has_value());
  ASSERT_FALSE(status.delta.has_value());
}

TEST_F(ConvergenceFinalFluxOrNTest, BadMaxIterations) {
  bart::convergence::FinalFluxOrN flux_tester(std::move(mock_multi_checker_ptr),
                                              system_fluxes);
  EXPECT_ANY_THROW(flux_tester.SetMaxIterations(0));
  EXPECT_ANY_THROW(flux_tester.SetMaxIterations(-1));
}

TEST_F(ConvergenceFinalFluxOrNTest, BadIterations) {
  bart::convergence::FinalFluxOrN flux_tester(std::move(mock_multi_checker_ptr),
                                              system_fluxes);
  EXPECT_ANY_THROW(flux_tester.SetIteration(0));
  EXPECT_ANY_THROW(flux_tester.SetIteration(-1));
}

/* Checks for proper handling of the MultiChecker returning a status of
 * complete, without reaching max iterations */
TEST_F(ConvergenceFinalFluxOrNTest, Convergence) {
  EXPECT_CALL(*mock_multi_checker_ptr, CheckIfConverged(_,_)).
      WillOnce(::testing::Return(true));

  bart::convergence::FinalFluxOrN flux_tester(std::move(mock_multi_checker_ptr),
                                              system_fluxes);
  flux_tester.SetMaxIterations(5).SetIteration(0);
  auto status = flux_tester.CheckFinalConvergence();
  EXPECT_TRUE(status.is_complete);
  EXPECT_EQ(status.iteration_number, 1);
  EXPECT_FALSE(status.failed_index.has_value());
  EXPECT_FALSE(status.delta.has_value());
}

/* Checks for proper handling of the Multichecker returning a status of
 * not complete, without reaching max iterations */
TEST_F(ConvergenceFinalFluxOrNTest, NoConvergence) {
  EXPECT_CALL(*mock_multi_checker_ptr, CheckIfConverged(_,_)).
      WillOnce(::testing::Return(false));
  EXPECT_CALL(*mock_multi_checker_ptr, failed_index()).
      WillOnce(::testing::Return(std::make_optional(2)));
  EXPECT_CALL(*mock_multi_checker_ptr, failed_delta()).
      WillOnce(::testing::Return(std::make_optional(1.3)));

  bart::convergence::FinalFluxOrN flux_tester(std::move(mock_multi_checker_ptr),
                                              system_fluxes);
  flux_tester.SetMaxIterations(5).SetIteration(0);
  auto status = flux_tester.CheckFinalConvergence();
  EXPECT_FALSE(status.is_complete);
  EXPECT_EQ(status.iteration_number, 1);
  EXPECT_EQ(status.failed_index.value(), 2);
  EXPECT_DOUBLE_EQ(status.delta.value(), 1.3);
}

/* Checks for proper handling of the Multichecker returning a status of
 * complete when max iterations is reached, even if not converged. */
TEST_F(ConvergenceFinalFluxOrNTest, MaxIterNoConvergence) {
  EXPECT_CALL(*nice_mock_multi_checker_ptr, CheckIfConverged(_,_)).
      WillOnce(::testing::Return(false));
  ON_CALL(*nice_mock_multi_checker_ptr, failed_index()).
      WillByDefault(::testing::Return(std::make_optional(2)));
  ON_CALL(*nice_mock_multi_checker_ptr, failed_delta()).
      WillByDefault(::testing::Return(std::make_optional(1.3)));

  bart::convergence::FinalFluxOrN flux_tester(
      std::move(nice_mock_multi_checker_ptr),
      system_fluxes);

  flux_tester.SetMaxIterations(5).SetIteration(4);
  auto status = flux_tester.CheckFinalConvergence();
  EXPECT_TRUE(status.is_complete);
}

/* Checks for proper handling of the Multichecker returning a status of
 * complete when max iterations is reached, if converged. */
TEST_F(ConvergenceFinalFluxOrNTest, MaxIterConvergence) {
  EXPECT_CALL(*mock_multi_checker_ptr, CheckIfConverged(_,_)).
      WillOnce(::testing::Return(true));

  bart::convergence::FinalFluxOrN flux_tester(std::move(mock_multi_checker_ptr),
                                              system_fluxes);

  flux_tester.SetMaxIterations(5).SetIteration(4);
  auto status = flux_tester.CheckFinalConvergence();
  EXPECT_TRUE(status.is_complete);
}

