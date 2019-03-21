#include "convergence/moments/multi_moment_checker_max.h"

#include <memory>

#include "convergence/moments/tests/single_moment_checker_mock.h"
#include "data/moment_types.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

using ::testing::_;
using ::testing::Return;
using ::testing::ByRef;

class MultiMomentCheckerMaxTest : public ::testing::Test {
 protected:
  using MultiMomentCheckerMax =
      bart::convergence::moments::MultiMomentCheckerMax;
  using SingleMomentCheckerMock =
      bart::convergence::moments::SingleMomentCheckerMock;

  bart::data::MomentsMap moments_map_one;
  bart::data::MomentsMap moments_map_two;

  std::unique_ptr<SingleMomentCheckerMock> checker_ptr;

  void SetUp() override;
};

void MultiMomentCheckerMaxTest::SetUp() {
  // Fill the two moments maps with all the harmonics for g = 5, l = 4
  int groups = 5;
  int max_l = 4;
  for (int group = 0; group < groups; ++group) {
    for (int l = 0; l <= max_l; ++l) {
      for (int m = -l; m <= max_l; ++m) {
        auto random_vector = btest::RandomVector(5, 0, 10);
        bart::data::MomentVector temp_moment(random_vector.cbegin(),
                                             random_vector.cend());
        moments_map_one[{group, l, m}] = temp_moment;
        moments_map_two[{group, l, m}] = temp_moment;
      }
    }
  }

  checker_ptr =
      std::make_unique<::testing::NiceMock<SingleMomentCheckerMock>>();

  // Single checker mock will by default return true about status of convergence
  ON_CALL(*checker_ptr, CheckIfConverged(_,_))
      .WillByDefault(Return(true));
}

TEST_F(MultiMomentCheckerMaxTest, Constructor) {
  MultiMomentCheckerMax test_checker(std::move(checker_ptr));
  EXPECT_EQ(checker_ptr, nullptr);
}

// -- MATCHING --

TEST_F(MultiMomentCheckerMaxTest, GoodMatch) {
  MultiMomentCheckerMax test_checker(std::move(checker_ptr));
  EXPECT_TRUE(test_checker.CheckIfConverged(moments_map_one, moments_map_two));
  EXPECT_TRUE(test_checker.is_converged());
  EXPECT_EQ(test_checker.failed_index(), std::nullopt);
  EXPECT_EQ(test_checker.delta(), std::nullopt);
}

TEST_F(MultiMomentCheckerMaxTest, BadMatch) {
  int failing_group = 2;
  // Set up mock to indicate that a specific group will fail
  bart::data::MomentVector failing_moment =
      moments_map_two[{failing_group, 0, 0}];

  EXPECT_CALL(*checker_ptr,
      CheckIfConverged(failing_moment, failing_moment))
      .WillOnce(Return(false));
  EXPECT_CALL(*checker_ptr, delta())
      .WillOnce(Return(0.123));

  MultiMomentCheckerMax test_checker(std::move(checker_ptr));

  EXPECT_TRUE(test_checker.CheckIfConverged(moments_map_one, moments_map_two));
  EXPECT_TRUE(test_checker.is_converged());
  EXPECT_EQ(test_checker.failed_index(), failing_group);
  EXPECT_EQ(test_checker.delta(), 0.123);
}

// -- ERRORS --

// Passing empty moments map should throw an error
TEST_F(MultiMomentCheckerMaxTest, EmptyMoments) {
  MultiMomentCheckerMax test_checker(std::move(checker_ptr));
  bart::data::MomentsMap empty_map;
  EXPECT_ANY_THROW(test_checker.CheckIfConverged(empty_map, empty_map));
}

// Passing a moment maps where one is missing a group should throw an error
TEST_F(MultiMomentCheckerMaxTest, WrongGroups) {
  MultiMomentCheckerMax test_checker(std::move(checker_ptr));
  auto iterator_to_delete = moments_map_one.find({2, 0, 0});
  moments_map_one.erase(iterator_to_delete);
  EXPECT_ANY_THROW(
      test_checker.CheckIfConverged(moments_map_one, moments_map_two));

}