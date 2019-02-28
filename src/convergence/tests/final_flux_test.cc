#include "convergence/final_flux.h"

#include <memory>
#include <optional>

#include "convergence/status.h"
#include "convergence/flux/tests/multi_checker_mock.h"
#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"

using ::testing::_;

class ConvergenceFinalFluxTest : public ::testing::Test {
 protected:
  bart::convergence::FinalFlux flux_tester;
  std::unique_ptr<bart::convergence::flux::MultiCheckerI> multi_checker_ptr;
  std::unique_ptr<bart::convergence::flux::MultiCheckerMock> mock_multi_checker;
  void SetUp();
  void MockToPointer();
};

void ConvergenceFinalFluxTest::SetUp() {
  mock_multi_checker =
      std::make_unique<bart::convergence::flux::MultiCheckerMock>();
}

void ConvergenceFinalFluxTest::MockToPointer() {
  multi_checker_ptr = std::move(mock_multi_checker);
}

TEST_F(ConvergenceFinalFluxTest, DefaultStatus) {
  bart::convergence::Status status;
  status = flux_tester.convergence_status();
  ASSERT_EQ(status.iteration_number, 0);
  ASSERT_EQ(status.max_iterations, 1);
  ASSERT_FALSE(status.is_complete);
  ASSERT_FALSE(status.failed_index.has_value());
  ASSERT_FALSE(status.delta.has_value());
}

