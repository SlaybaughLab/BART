#include "../group_flux_checker_sequential.h"

#include <memory>

#include <deal.II/lac/petsc_parallel_vector.h>

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

// Tests where there are no mock calls
class GroupFluxCheckerSeqTestEmptyMock : public GroupFluxCheckerSequentialTest {
 protected:
  bart::convergence::GroupFluxCheckerSequential sequential_tester;
  void SetUp() override {
    MocksToPointers();
    sequential_tester.ProvideChecker(tester_ptr);
  }
};

TEST_F(GroupFluxCheckerSeqTestEmptyMock, Constructor) {
  EXPECT_EQ(tester_ptr, nullptr);
}

TEST_F(GroupFluxCheckerSeqTestEmptyMock, DifferentGroupSizes) {
  EXPECT_TRUE(true);
}




