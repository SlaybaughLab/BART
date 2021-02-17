#include "calculator/residual/isotropic_residual.hpp"

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include "calculator/residual/tests/vector_difference_mock.hpp"
#include "system/moments/tests/spherical_harmonic_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_assertions.hpp"
#include "test_helpers/test_helper_functions.h"

#include "calculator/residual/tests/isotropic_residual_mock.hpp"

namespace  {

using namespace bart;
using ::testing::NiceMock, ::testing::Return, ::testing::DoDefault, ::testing::ReturnRef, ::testing::Ref;

class CalculatorResidualIsotropicResidualTest : public ::testing::Test {
 public:
  using IsotropicResidualCalculator = calculator::residual::IsotropicResidual;
  using ResidualCalculatorMock = NiceMock<calculator::residual::VectorDifferenceMock>;
  using FullMatrix = dealii::FullMatrix<double>;
  using MomentsMock = system::moments::SphericalHarmonicMock;
  using Vector = dealii::Vector<double>;

  // Test object
  std::unique_ptr<IsotropicResidualCalculator> test_calculator_;

  // Supporting test objects
  std::shared_ptr<MomentsMock> half_step_moments_{ std::make_shared<MomentsMock>() };
  std::shared_ptr<MomentsMock> previous_step_moments_{ std::make_shared<MomentsMock>() };
  std::map<int, Vector> half_step_moment_vectors_, previous_step_moment_vectors_, residual_vectors_;

  // Mock observation pointers
  ResidualCalculatorMock* residual_calculator_mock_obs_ptr_;

  // Test parameters
  FullMatrix sigma_s_;
  Vector expected_isotropic_residual_;
  const int total_groups_{ test_helpers::RandomInt(5, 10) };
  const int test_group_{ test_helpers::RandomInt(0, total_groups_ - 2)};
  const int vector_size_{ test_helpers::RandomInt(10, 20) };
  auto SetUp() -> void override;
};

auto CalculatorResidualIsotropicResidualTest::SetUp() -> void {
  auto residual_calculator_mock_ptr = std::make_unique<ResidualCalculatorMock>();
  residual_calculator_mock_obs_ptr_ = residual_calculator_mock_ptr.get();
  test_calculator_ = std::make_unique<IsotropicResidualCalculator>(std::move(residual_calculator_mock_ptr));

  auto get_vector = [](const int size, const double value_to_set) {
    Vector return_vector(size);
    for (auto& val : return_vector)
      val = value_to_set;
    return return_vector;
  };

  sigma_s_.reinit(total_groups_, total_groups_);
  for (int g = 0; g < total_groups_; ++g) {
    for (int g_in = g + 1; g_in < total_groups_; ++g_in) {
      sigma_s_(g, g_in) = g + g_in + 1;
    }
  }

  for (int group = 0; group < total_groups_; ++group) {
    std::array<int, 3> index{group, 0, 0};
    half_step_moment_vectors_[group] = get_vector(vector_size_, 1.5*group);
    previous_step_moment_vectors_[group] = get_vector(vector_size_, 0.5*group);

    ON_CALL(*half_step_moments_, GetMoment(index)).WillByDefault(ReturnRef(half_step_moment_vectors_.at(group)));
    ON_CALL(*previous_step_moments_, GetMoment(index)).WillByDefault(ReturnRef(previous_step_moment_vectors_.at(group)));
  }
  ON_CALL(*half_step_moments_, total_groups()).WillByDefault(Return(total_groups_));
  ON_CALL(*previous_step_moments_, total_groups()).WillByDefault(Return(total_groups_));

  expected_isotropic_residual_.reinit(vector_size_);

  double isotropic_residual_value{ 0 };
  for (int group = test_group_ + 1; group < total_groups_; ++group) {
    const double residual_value = (1 + test_group_ + group) * group;
    isotropic_residual_value += residual_value;

    residual_vectors_[group] = get_vector(vector_size_, residual_value);
    ON_CALL(*residual_calculator_mock_obs_ptr_, CalculateResidual(Ref(half_step_moment_vectors_.at(group)),
                                                                  Ref(previous_step_moment_vectors_.at(group)),
                                                                  sigma_s_(test_group_, group)))
        .WillByDefault(Return(residual_vectors_.at(group)));
  }
  for (auto& val : expected_isotropic_residual_)
    val = isotropic_residual_value;
}

// test calculator should calculate the correct isotropic residual based on test values
TEST_F(CalculatorResidualIsotropicResidualTest, CalculateIsotropicResidual) {
  EXPECT_CALL(*previous_step_moments_, total_groups()).WillOnce(DoDefault());
  EXPECT_CALL(*half_step_moments_, total_groups()).WillOnce(DoDefault());

  for (int group_in = test_group_ + 1; group_in < total_groups_; ++group_in) {
    std::array<int, 3> index{ group_in, 0, 0 };
    EXPECT_CALL(*previous_step_moments_, GetMoment(index)).WillOnce(DoDefault());
    EXPECT_CALL(*half_step_moments_, GetMoment(index)).WillOnce(DoDefault());
    EXPECT_CALL(*residual_calculator_mock_obs_ptr_, CalculateResidual(
        Ref(half_step_moment_vectors_.at(group_in)), Ref(previous_step_moment_vectors_.at(group_in)),
        sigma_s_(test_group_, group_in)))
        .WillOnce(DoDefault());
  }

  const auto isotropic_residual = test_calculator_->CalculateIsotropicResidual(half_step_moments_.get(),
                                                                               previous_step_moments_.get(),
                                                                               test_group_, sigma_s_);
  EXPECT_TRUE(test_helpers::AreEqual(isotropic_residual, expected_isotropic_residual_));
}

// If two moment objects return different total group values, throw
TEST_F(CalculatorResidualIsotropicResidualTest, BadTotalGroupsInMoments) {
  EXPECT_CALL(*previous_step_moments_, total_groups()).WillOnce(Return(total_groups_ + 1));
  EXPECT_CALL(*half_step_moments_, total_groups()).WillOnce(DoDefault());
  EXPECT_ANY_THROW({
    [[maybe_unused]] const auto isotropic_residual =
        test_calculator_->CalculateIsotropicResidual(half_step_moments_.get(), previous_step_moments_.get(),
                                                     test_group_, sigma_s_);
  });
  EXPECT_CALL(*previous_step_moments_, total_groups()).WillOnce(DoDefault());
  EXPECT_CALL(*half_step_moments_, total_groups()).WillOnce(Return(total_groups_ + 1));
  EXPECT_ANY_THROW({
    [[maybe_unused]] const auto isotropic_residual =
        test_calculator_->CalculateIsotropicResidual(half_step_moments_.get(), previous_step_moments_.get(),
                                                    test_group_, sigma_s_);
  });
}

// If specified group number is negative, throw
TEST_F(CalculatorResidualIsotropicResidualTest, NegativeGroup) {
  const std::vector<int> bad_groups{ -1, -2, -100 };
  for (const auto bad_group : bad_groups){
    EXPECT_ANY_THROW({
      [[maybe_unused]] const auto isotropic_residual =
          test_calculator_->CalculateIsotropicResidual(half_step_moments_.get(), previous_step_moments_.get(),
                                                       bad_group, sigma_s_);
    });
  }
}

// If specified group exceeds total groups, throw
TEST_F(CalculatorResidualIsotropicResidualTest, GroupExceedsTotalGroups) {
  const std::vector<int> bad_groups{ total_groups_, total_groups_ + 1, total_groups_ + 100 };
  for (const auto bad_group : bad_groups){
    EXPECT_CALL(*previous_step_moments_, total_groups()).WillOnce(DoDefault());
    EXPECT_CALL(*half_step_moments_, total_groups()).WillOnce(DoDefault());
    EXPECT_ANY_THROW({
      [[maybe_unused]] const auto isotropic_residual =
          test_calculator_->CalculateIsotropicResidual(half_step_moments_.get(), previous_step_moments_.get(),
                                                       bad_group, sigma_s_);
                     });
  }
}



} // namespace
