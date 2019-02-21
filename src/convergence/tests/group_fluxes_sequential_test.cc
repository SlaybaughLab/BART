#include "../group_fluxes_sequential.h"

#include <memory>

#include "../../test_helpers/gmock_wrapper.h"
#include "flux_checker_mock.h"

class GroupFluxesSequentialTest : public ::testing::Test {
 protected:
  std::unique_ptr<bart::convergence::FluxCheckerI> tester_ptr;
  std::unique_ptr<bart::convergence::FluxCheckerMock> tester_mock;

  void SetUp() override;
  void MocksToPointers() {
    tester_ptr = std::move(tester_mock);
  };
};

void GroupFluxesSequentialTest::SetUp() {
  tester_mock = std::make_unique<bart::convergence::FluxCheckerMock>();
}

TEST_F(GroupFluxesSequentialTest, Constructor) {
  MocksToPointers();
  bart::convergence::GroupFluxesSequential sequential_tester(tester_ptr);
  EXPECT_EQ(tester_ptr, nullptr);
}



