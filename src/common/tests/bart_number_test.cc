#include "../numbers.h"

#include "gtest/gtest.h"

class NumbersTest : public ::testing::Test {
 protected:
  void SetUp ();

  void ConstNumbersTest ();

  double pi_;
};

void NumbersTest::SetUp() {
  pi_ = dealii::numbers::PI;
}

TEST_F(NumbersTest, ConstNumTest) {
  EXPECT_FLOAT_EQ(pi_, bconst::kPi);
  EXPECT_FLOAT_EQ(2.*pi_, bconst::kTwoPi);
  EXPECT_FLOAT_EQ(4.*pi_, bconst::kFourPi);
  EXPECT_FLOAT_EQ(1./pi_, bconst::kInvPi);
  EXPECT_FLOAT_EQ(0.5/pi_, bconst::kInvTwoPi);
  EXPECT_FLOAT_EQ(0.25/pi_, bconst::kInvFourPi);
  EXPECT_FLOAT_EQ(1.0e-15, bconst::kSmall);
}
