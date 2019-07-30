#include "iteration/outer/outer_power_iteration.h"

#include <memory>

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class IterationOuterPowerIterationTest : public ::testing::Test {
 protected:
  static constexpr int dim = DimensionWrapper::value;

  using OuterPowerIteration = iteration::outer::OuterPowerIteration;
  std::unique_ptr<OuterPowerIteration> test_iterator;

  void SetUp() override;
};

template <typename DimensionWrapper>
void IterationOuterPowerIterationTest<DimensionWrapper>::SetUp() {
  test_iterator = std::make_unique<OuterPowerIteration>();
}


TYPED_TEST_CASE(IterationOuterPowerIterationTest, bart::testing::AllDimensions);

TYPED_TEST(IterationOuterPowerIterationTest, Constructor) {
  EXPECT_NE(this->test_iterator, nullptr);
}

} // namespace