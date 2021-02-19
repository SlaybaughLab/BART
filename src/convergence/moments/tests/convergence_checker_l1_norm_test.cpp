#include "convergence/moments/convergence_checker_l1_norm.hpp"

#include "system/moments/spherical_harmonic_types.h"
#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

class SingleMomentCheckerL1NormTest : public ::testing::Test {
 public:
  using MomentVector = bart::system::moments::MomentVector;
  static constexpr double max_delta{ 1e-6 };
  static constexpr MomentVector::size_type moment_size{ 5 };
  bart::convergence::moments::ConvergenceCheckerL1Norm checker{max_delta };
  bart::system::moments::MomentVector moment_one, moment_two;
  auto SetUp() -> void override;
};

auto SingleMomentCheckerL1NormTest::SetUp() -> void {
  auto set_up_moment = [=](MomentVector& moment) {
    moment.reinit(moment_size);
    const auto random_vector{ test_helpers::RandomVector(5, 0, 2) };
    for (MomentVector::size_type i = 0; i < moment_size; ++i)
      moment[i] = random_vector.at(i);
  };
  set_up_moment(moment_one);
  moment_two = moment_one;
}

// Max Delta getter should return the correct value
TEST_F(SingleMomentCheckerL1NormTest, Getters) {
  EXPECT_EQ(checker.max_delta(), this->max_delta);
}

// Max delta setter should set value correctly
TEST_F(SingleMomentCheckerL1NormTest, MaxDeltaSetter) {
  const double max_delta_to_set{ test_helpers::RandomDouble(1e-16, 1e-8) };
  EXPECT_EQ(checker.max_delta(), this->max_delta);
  checker.SetMaxDelta(max_delta_to_set);
  EXPECT_EQ(checker.max_delta(), max_delta_to_set);
}

// Setter should throw if provided a negative value
TEST_F(SingleMomentCheckerL1NormTest, SetterBadValue) {
  const double max_delta_to_set{ test_helpers::RandomDouble(-1e-8, -1e-16) };
  EXPECT_EQ(checker.max_delta(), this->max_delta);
  EXPECT_ANY_THROW({checker.SetMaxDelta(max_delta_to_set);});
}

// The same vector should return true
TEST_F(SingleMomentCheckerL1NormTest, SameVector) {
  EXPECT_TRUE(checker.IsConverged(moment_one, moment_one));
  EXPECT_TRUE(checker.is_converged());
}

// Being slightly less than one max delta away should return true
TEST_F(SingleMomentCheckerL1NormTest, OneThresholdAway) {
  double to_add = moment_one.l1_norm() * 0.99 * checker.max_delta();
  moment_two(2) += to_add;

  EXPECT_TRUE(checker.IsConverged(moment_one, moment_two));
  EXPECT_TRUE(checker.IsConverged(moment_two, moment_one));
  EXPECT_TRUE(checker.is_converged());
  EXPECT_NEAR(0.99 * checker.max_delta(), checker.delta().value(), 1e-6);
}

// Being greater than max delta away should return false
TEST_F(SingleMomentCheckerL1NormTest, TwoThresholdAway) {
  const double to_add = moment_one.l1_norm() * 2 * checker.max_delta();
  moment_two(2) += to_add;

  EXPECT_FALSE(checker.IsConverged(moment_one, moment_two));
  EXPECT_FALSE(checker.IsConverged(moment_two, moment_one));
  EXPECT_FALSE(checker.is_converged());
  EXPECT_NEAR(2 * checker.max_delta(), checker.delta().value(), 1e-6);
}

// Setting new max delta should still work for convergence
TEST_F(SingleMomentCheckerL1NormTest, SetMaxDelta) {
  const double to_set{ 1e-5 };
  checker.SetMaxDelta(to_set);
  EXPECT_EQ(checker.max_delta(), to_set);

  const double to_add = moment_one.l1_norm() * 0.99 * to_set;
  moment_two(2) += to_add;

  EXPECT_TRUE(checker.IsConverged(moment_one, moment_two));
  EXPECT_TRUE(checker.IsConverged(moment_two, moment_one));
}

} // namespace