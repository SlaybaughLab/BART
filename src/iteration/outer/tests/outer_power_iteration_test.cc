#include "iteration/outer/outer_power_iteration.h"

#include <memory>

#include "iteration/updater/tests/source_updater_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "system/system.h"

namespace  {

using namespace bart;

using ::testing::Ref;

template <typename DimensionWrapper>
class IterationOuterPowerIterationTest : public ::testing::Test {
 protected:
  static constexpr int dim = DimensionWrapper::value;
  using OuterPowerIteration = iteration::outer::OuterPowerIteration;
  using SourceUpdater = iteration::updater::SourceUpdaterMock;

  std::unique_ptr<OuterPowerIteration> test_iterator;

  // Dependencies
  std::shared_ptr<SourceUpdater> source_updater_ptr_;

  // Supporting objects
  system::System test_system;

  // Observation pointers

  // Test parameters
  const int total_groups = 2;
  const int total_angles = 3;
  const int iterations_ = 4;

  void SetUp() override;
};

template <typename DimensionWrapper>
void IterationOuterPowerIterationTest<DimensionWrapper>::SetUp() {
  source_updater_ptr_ = std::make_shared<SourceUpdater>();

  // Set up system
  test_system.total_angles = total_angles;
  test_system.total_groups = total_groups;

  // Construct test object
  test_iterator = std::make_unique<OuterPowerIteration>(
      source_updater_ptr_
      );
}


TYPED_TEST_CASE(IterationOuterPowerIterationTest, bart::testing::AllDimensions);

TYPED_TEST(IterationOuterPowerIterationTest, Constructor) {
  EXPECT_NE(this->test_iterator, nullptr);
  EXPECT_NE(this->test_iterator->source_updater_ptr(), nullptr);
  EXPECT_EQ(this->source_updater_ptr_.use_count(), 2);
}

TYPED_TEST(IterationOuterPowerIterationTest, ConstructorErrors) {
  EXPECT_ANY_THROW({
    iteration::outer::OuterPowerIteration test_iterator(nullptr);
  });
}

TYPED_TEST(IterationOuterPowerIterationTest, IterateToConvergenceTest) {
  for (int group = 0; group < this->total_groups; ++group) {
    for (int angle = 0; angle < this->total_angles; ++angle) {
      EXPECT_CALL(*this->source_updater_ptr_, UpdateFissionSource(
          Ref(this->test_system),group, angle))
          .Times(this->iterations_);
    }
  }

  this->test_iterator->IterateToConvergence(this->test_system);
}



} // namespace