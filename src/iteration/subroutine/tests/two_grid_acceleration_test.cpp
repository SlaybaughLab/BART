#include "iteration/subroutine/two_grid_acceleration.hpp"

#include <deal.II/lac/vector.h>

#include "acceleration/two_grid/tests/flux_corrector_mock.hpp"
#include "calculator/residual/tests/isotropic_residual_mock.hpp"
#include "data/cross_sections/tests/cross_sections_mock.hpp"
#include "framework/tests/framework_mock.hpp"
#include "system/moments/tests/spherical_harmonic_mock.h"
#include "system/system.hpp"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class IterationSubroutineTwoGridAccelerationTest : public ::testing::Test {
 public:
  using CrossSectionsMock = data::cross_sections::CrossSectionsMock;
  using RHSVector = dealii::Vector<double>;
  using FluxCorrectorMock = acceleration::two_grid::FluxCorrectorMock;
  using FrameworkMock = framework::FrameworkMock;
  using IsotropicResidualCalculatorMock = calculator::residual::IsotropicResidualMock;
  using TwoGridSubroutine = iteration::subroutine::TwoGridAcceleration;
  using System = system::System;

  // Dependency pointers
  std::shared_ptr<RHSVector> rhs_vector_ptr_;
  std::shared_ptr<CrossSectionsMock> cross_sections_ptr_ { std::make_shared<CrossSectionsMock>() };
  FluxCorrectorMock* flux_corrector_mock_obs_ptr_;
  FrameworkMock* framework_mock_obs_ptr_;
  IsotropicResidualCalculatorMock* isotropic_residual_calculator_mock_obs_ptr_;

  std::unique_ptr<TwoGridSubroutine> test_subroutine_;
  System test_system_;

  static constexpr int vector_size_{ 10 };
  static constexpr int total_groups_{ 3 };
  auto SetUp() -> void override;
};

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
                                                         cross_sections_ptr_,
                                                         rhs_vector_ptr_);

  test_system_.total_groups = this->total_groups_;
}

// Getters for depdendencies should return correct pointers to dependencies
TEST_F(IterationSubroutineTwoGridAccelerationTest, DependencyGetters) {
  EXPECT_EQ(test_subroutine_->framework_ptr(), framework_mock_obs_ptr_);
  EXPECT_EQ(test_subroutine_->flux_corrector_ptr(), flux_corrector_mock_obs_ptr_);
  EXPECT_EQ(test_subroutine_->residual_calculator_ptr(), isotropic_residual_calculator_mock_obs_ptr_);
  EXPECT_EQ(test_subroutine_->cross_sections_ptr(), cross_sections_ptr_);
  EXPECT_EQ(test_subroutine_->isotropic_residual_ptr(), rhs_vector_ptr_);
}

// Passing a null dependency should throw an error
TEST_F(IterationSubroutineTwoGridAccelerationTest, NullDependenciesThrow) {
  constexpr int n_dependencies{ 5 };
  for (int i = 0; i < n_dependencies; ++i) {
    EXPECT_ANY_THROW({
      TwoGridSubroutine(i == 0 ? nullptr : std::make_unique<FluxCorrectorMock>(),
                        i == 1 ? nullptr : std::make_unique<FrameworkMock>(),
                        i == 2 ? nullptr : std::make_unique<IsotropicResidualCalculatorMock>(),
                        i == 3 ? nullptr : cross_sections_ptr_,
                        i == 4 ? nullptr : rhs_vector_ptr_);
    });
  }
}

TEST_F(IterationSubroutineTwoGridAccelerationTest, Execute) {

}

} // namespace
