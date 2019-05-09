#include "data/system/left_hand_side_fixed.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class SystemLeftHandSideFixedTest : public ::testing::Test {
  data::system::LeftHandSideFixed test_lhs_;
};

TEST_F(SystemLeftHandSideFixedTest, Dummy) {
  EXPECT_TRUE(true);
}

} // namespace