#include "results/output_dealii_vtu.h"

#include "test_helpers/gmock_wrapper.h"

namespace {

template <typename DimensionWrapper>
class OutputDealiiVtuTest : public ::testing::Test {

};

TYPED_TEST_CASE(OutputDealiiVtuTest, bart::testing::AllDimensions);

TYPED_TEST(OutputDealiiVtuTest, Constructor) {
  EXPECT_TRUE(true);
}

}