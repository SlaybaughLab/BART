#include "../group_flux_checker_sequential.h"

#include <memory>

#include "../../test_helpers/gmock_wrapper.h"
#include "flux_checker_mock.h"

class GroupFluxCheckerSequentialTest : public ::testing::Test {
 protected:
  std::unique_ptr<bart::convergence::FluxCheckerI> tester_ptr;
  std::unique_ptr<bart::convergence::FluxCheckerMock> tester_mock;

  void SetUp() override;
  void MocksToPointers() {
    tester_ptr = std::move(tester_mock);
  };
};

void GroupFluxCheckerSequentialTest::SetUp() {
  tester_mock = std::make_unique<bart::convergence::FluxCheckerMock>();
}

TEST_F(GroupFluxCheckerSequentialTest, Constructor) {
  MocksToPointers();
  bart::convergence::GroupFluxCheckerSequential sequential_tester(tester_ptr);
  EXPECT_EQ(tester_ptr, nullptr);
}



