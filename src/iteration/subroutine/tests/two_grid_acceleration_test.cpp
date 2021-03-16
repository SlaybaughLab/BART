#include "iteration/subroutine/two_grid_acceleration.hpp"

#include <deal.II/lac/vector.h>

#include "acceleration/two_grid/tests/flux_corrector_mock.hpp"
#include "calculator/residual/tests/domain_isotropic_residual_mock.hpp"
#include "framework/tests/framework_mock.hpp"
#include "system/moments/tests/spherical_harmonic_mock.h"
#include "system/system.hpp"
#include "test_helpers/test_helper_functions.h"
#include "test_helpers/test_assertions.hpp"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;
using ::testing::NiceMock, ::testing::DoDefault, ::testing::Return, ::testing::ReturnRef, ::testing::_;
using ::testing::Ref, ::testing::AtLeast;

class IterationSubroutineTwoGridAccelerationTest : public ::testing::Test {
 public:
  using RHSVector = dealii::Vector<double>;
  using FluxCorrectorMock = acceleration::two_grid::FluxCorrectorMock;
  using FrameworkMock = framework::FrameworkMock;
  using IsotropicResidualCalculatorMock = calculator::residual::DomainIsotropicResidualMock;
  using Moments = system::moments::SphericalHarmonicMock;
  using TwoGridSubroutine = iteration::subroutine::TwoGridAcceleration;
  using System = system::System;

  // Dependency pointers
  std::shared_ptr<RHSVector> rhs_vector_ptr_;
  FluxCorrectorMock* flux_corrector_mock_obs_ptr_;
  FrameworkMock* framework_mock_obs_ptr_;
  IsotropicResidualCalculatorMock* isotropic_residual_calculator_mock_obs_ptr_;
  Moments* current_moments_obs_ptr_;
  Moments* previous_moments_obs_ptr_;
  Moments* framework_system_current_moments_obs_ptr_;

  std::unique_ptr<TwoGridSubroutine> test_subroutine_;
  System test_system_;
  std::shared_ptr<System> framework_system_{ std::make_shared<System>() };

  dealii::Vector<double> total_isotropic_residual_;

  static constexpr int vector_size_{ 10 };
  static constexpr int total_groups_{ 3 };
  auto SetUp() -> void override;
};

auto DealiiVector(std::vector<double> to_convert) {
  dealii::Vector<double> return_vector(to_convert.size());
  for (std::vector<double>::size_type i = 0; i < to_convert.size(); ++i)
    return_vector(i) = to_convert.at(i);
  return return_vector;
}

auto IterationSubroutineTwoGridAccelerationTest::SetUp() -> void {
  rhs_vector_ptr_ = std::make_shared<RHSVector>(vector_size_);

  auto flux_corrector_ptr = std::make_unique<FluxCorrectorMock>();
  flux_corrector_mock_obs_ptr_ = flux_corrector_ptr.get();
  auto isotropic_residual_calculator_ptr = std::make_unique<IsotropicResidualCalculatorMock>();
  isotropic_residual_calculator_mock_obs_ptr_ = isotropic_residual_calculator_ptr.get();
  auto framework_ptr = std::make_unique<FrameworkMock>();
  framework_mock_obs_ptr_ = framework_ptr.get();

  test_subroutine_ = std::make_unique<TwoGridSubroutine>(std::move(flux_corrector_ptr),
                                                         std::move(framework_ptr),
                                                         std::move(isotropic_residual_calculator_ptr),
                                                         rhs_vector_ptr_);

  test_system_.total_groups = this->total_groups_;
  test_system_.current_moments = std::make_shared<Moments>();
  test_system_.previous_moments = std::make_unique<Moments>();
  current_moments_obs_ptr_ = dynamic_cast<Moments*>(test_system_.current_moments.get());
  previous_moments_obs_ptr_ = dynamic_cast<Moments*>(test_system_.previous_moments.get());

  framework_system_->total_groups = 1;
  framework_system_->current_moments = std::make_shared<Moments>();
  framework_system_current_moments_obs_ptr_ = dynamic_cast<Moments*>(framework_system_->current_moments.get());

  total_isotropic_residual_ = DealiiVector(test_helpers::RandomVector(vector_size_, 0, 100));
  ON_CALL(*isotropic_residual_calculator_mock_obs_ptr_, CalculateDomainResidual(_,_))
      .WillByDefault(Return(total_isotropic_residual_));
}

// Getters for depdendencies should return correct pointers to dependencies
TEST_F(IterationSubroutineTwoGridAccelerationTest, DependencyGetters) {
  EXPECT_EQ(test_subroutine_->framework_ptr(), framework_mock_obs_ptr_);
  EXPECT_EQ(test_subroutine_->flux_corrector_ptr(), flux_corrector_mock_obs_ptr_);
  EXPECT_EQ(test_subroutine_->residual_calculator_ptr(), isotropic_residual_calculator_mock_obs_ptr_);
  EXPECT_EQ(test_subroutine_->isotropic_residual_ptr(), rhs_vector_ptr_);
}

// Passing a null dependency should throw an error
TEST_F(IterationSubroutineTwoGridAccelerationTest, NullDependenciesThrow) {
  constexpr int n_dependencies{ 4 };
  for (int i = 0; i < n_dependencies; ++i) {
    EXPECT_ANY_THROW({
      TwoGridSubroutine(i == 0 ? nullptr : std::make_unique<FluxCorrectorMock>(),
                        i == 1 ? nullptr : std::make_unique<FrameworkMock>(),
                        i == 2 ? nullptr : std::make_unique<IsotropicResidualCalculatorMock>(),
                        i == 3 ? nullptr : rhs_vector_ptr_);
    });
  }
}

TEST_F(IterationSubroutineTwoGridAccelerationTest, Execute) {
  EXPECT_CALL(*this->isotropic_residual_calculator_mock_obs_ptr_, CalculateDomainResidual(current_moments_obs_ptr_,
                                                                                          previous_moments_obs_ptr_))
      .WillOnce(DoDefault());
  EXPECT_CALL(*this->framework_mock_obs_ptr_, SolveSystem());
  EXPECT_CALL(*this->framework_mock_obs_ptr_, system())
      .Times(AtLeast(1))
      .WillRepeatedly(Return(this->framework_system_.get()));
  dealii::Vector<double> error_vector(DealiiVector(test_helpers::RandomVector(this->vector_size_, 0, 100)));

  EXPECT_CALL(*this->framework_system_current_moments_obs_ptr_, GetMoment(std::array<int, 3>{0, 0, 0}))
      .WillOnce(ReturnRef(error_vector));
  for (int group = 0; group < this->total_groups_; ++group) {
    dealii::Vector<double> flux_to_correct(DealiiVector(test_helpers::RandomVector(this->vector_size_, 0, 100)));
    EXPECT_CALL(*this->current_moments_obs_ptr_, GetMoment(std::array<int, 3>{group, 0, 0}))
        .WillOnce(ReturnRef(flux_to_correct));
    EXPECT_CALL(*this->flux_corrector_mock_obs_ptr_, CorrectFlux(Ref(flux_to_correct),
                                                                 Ref(error_vector), group));
  }

  this->test_subroutine_->Execute(this->test_system_);

  EXPECT_TRUE(test_helpers::AreEqual(this->total_isotropic_residual_, *this->rhs_vector_ptr_));
}

} // namespace
