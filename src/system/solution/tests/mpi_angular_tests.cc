#include "system/solution/mpi_angular.h"

#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

class SolutionMPIAngularTests : public ::testing::Test {
 protected:
  system::solution::MPIAngular test_solution;
};

TEST_F(SolutionMPIAngularTests, Dummy) {
  EXPECT_TRUE(true);
}

} // namespace