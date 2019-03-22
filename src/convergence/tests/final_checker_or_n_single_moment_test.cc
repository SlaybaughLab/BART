#include "convergence/final_checker_or_n.h"

#include <memory>

#include <gtest/gtest.h>

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
} // namespace




