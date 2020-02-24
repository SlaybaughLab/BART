#include "convergence/moments/multi_moment_checker_max.h"

#include <memory>
#include <optional>

#include "convergence/moments/tests/single_moment_checker_mock.h"
#include "system/moments/spherical_harmonic_types.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

using namespace bart;

using ::testing::_, ::testing::Return, ::testing::ByRef, ::testing::Ne;
using ::testing::Sequence, ::testing::Expectation;

class MultiMomentCheckerMaxTest : public ::testing::Test {
 protected:
  using MultiMomentCheckerMax =
      bart::convergence::moments::MultiMomentCheckerMax;
  using SingleMomentCheckerMock =
      bart::convergence::moments::SingleMomentCheckerMock;

  bart::system::moments::MomentsMap moments_map_one;
  bart::system::moments::MomentsMap moments_map_two;

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
        auto random_vector = test_helpers::RandomVector(5, 0, 10);
        bart::system::moments::MomentVector temp_moment(random_vector.cbegin(),
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

/* Bad match, and with multiple failing groups, only the one with the highest
 * error will be recorded.
 */
TEST_F(MultiMomentCheckerMaxTest, BadMatch) {

  int failing_group = 2;

  bart::system::moments::MomentVector failing_moment =
      moments_map_two[{failing_group, 0, 0}];

  EXPECT_CALL(*checker_ptr,
              CheckIfConverged(Ne(failing_moment), Ne(failing_moment)))
              .WillRepeatedly(Return(false));
  EXPECT_CALL(*checker_ptr, delta())
      .WillRepeatedly(Return(std::make_optional<double>(0.1)));

  Expectation max_group =
      EXPECT_CALL(*checker_ptr,
          CheckIfConverged(failing_moment, failing_moment))
      .WillOnce(Return(false));
  Expectation max_delta =
      EXPECT_CALL(*checker_ptr, delta())
      .After(max_group)
      .WillOnce(Return(std::make_optional<double>(0.123)));

  EXPECT_CALL(*checker_ptr, delta())
      .After(max_delta)
      .WillRepeatedly(Return(std::make_optional<double>(0.1)));

  MultiMomentCheckerMax test_checker(std::move(checker_ptr));

  EXPECT_FALSE(test_checker.CheckIfConverged(moments_map_one, moments_map_two));
  EXPECT_FALSE(test_checker.is_converged());
  EXPECT_EQ(test_checker.failed_index().value_or(-1), failing_group);
  EXPECT_EQ(test_checker.delta().value_or(-1), 0.123);
}

// Good match after bad should clear delta and indices
TEST_F(MultiMomentCheckerMaxTest, ConvergeAfterBad) {

  int failing_group = 2;

  bart::system::moments::MomentVector failing_moment =
      moments_map_two[{failing_group, 0, 0}];

  EXPECT_CALL(*checker_ptr,
              CheckIfConverged(Ne(failing_moment), Ne(failing_moment)))
      .WillRepeatedly(Return(false));
  EXPECT_CALL(*checker_ptr, delta())
      .WillRepeatedly(Return(std::make_optional<double>(0.1)));

  Expectation max_group =
      EXPECT_CALL(*checker_ptr,
                  CheckIfConverged(failing_moment, failing_moment))
          .WillOnce(Return(false));
  EXPECT_CALL(*checker_ptr, delta())
      .After(max_group)
      .WillOnce(Return(std::make_optional<double>(0.123)));

  EXPECT_CALL(*checker_ptr, CheckIfConverged(_,_))
      .After(max_group)
      .WillRepeatedly(Return(true));

  MultiMomentCheckerMax test_checker(std::move(checker_ptr));

  EXPECT_FALSE(test_checker.CheckIfConverged(moments_map_one, moments_map_two));
  EXPECT_FALSE(test_checker.is_converged());
  EXPECT_EQ(test_checker.failed_index().value_or(-1), failing_group);
  EXPECT_EQ(test_checker.delta().value_or(-1), 0.123);

  EXPECT_TRUE(test_checker.CheckIfConverged(moments_map_one, moments_map_two));
  EXPECT_TRUE(test_checker.is_converged());
  EXPECT_EQ(test_checker.failed_index(), std::nullopt);
  EXPECT_EQ(test_checker.delta(), std::nullopt);
}

// -- ERRORS --

// Passing empty moments map should throw an error
TEST_F(MultiMomentCheckerMaxTest, EmptyMoments) {
  MultiMomentCheckerMax test_checker(std::move(checker_ptr));
  bart::system::moments::MomentsMap empty_map;
  EXPECT_ANY_THROW(test_checker.CheckIfConverged(empty_map, empty_map));
}

// Passing a moment maps where current is a different length should throw
TEST_F(MultiMomentCheckerMaxTest, WrongLength) {
  MultiMomentCheckerMax test_checker(std::move(checker_ptr));
  auto iterator_to_delete = moments_map_one.find({2, 1, 1});
  moments_map_one.erase(iterator_to_delete);
  EXPECT_ANY_THROW(
      test_checker.CheckIfConverged(moments_map_one, moments_map_two));
}

// Passing a moment maps where one is missing a group should throw an error
TEST_F(MultiMomentCheckerMaxTest, WrongGroup) {
  MultiMomentCheckerMax test_checker(std::move(checker_ptr));
  auto entry_to_change = moments_map_one.extract({2, 0, 0});
  entry_to_change.key() = {2, 10, 10};
  moments_map_one.insert(std::move(entry_to_change));
  EXPECT_ANY_THROW(
      test_checker.CheckIfConverged(moments_map_one, moments_map_two));
}

} // namespace