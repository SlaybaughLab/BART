#include "iteration/subroutine/get_scalar_flux_from_framework.hpp"

#include "iteration/subroutine/tests/subroutine_mock.hpp"
<<<<<<< HEAD
#include "test_helpers/gmock_wrapper.h"
=======
#include "system/system.hpp"
#include "system/moments/tests/spherical_harmonic_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"
>>>>>>> 6b54f712... Added implementation of Execute to GetScalarFluxFromFramework.

namespace  {

using namespace bart;
using ::testing::Return, ::testing::ReturnRef;

class GetScalarFluxFromFrameworkSubroutineTest : public ::testing::Test {
 public:
  using GetScalarFluxFromFramework = bart::iteration::subroutine::GetScalarFluxFromFramework;
  using Framework = framework::FrameworkMock;
  using SphericalHarmonicMoments = system::moments::SphericalHarmonicMock;
  using System = system::System;
  using Vector = dealii::Vector<double>;

  // Test objects
  System test_system_;
  std::shared_ptr<System> framework_system_ptr_{ nullptr };
  std::unique_ptr<GetScalarFluxFromFramework> test_subroutine_{ nullptr };
  std::map<int, Vector> framework_scalar_fluxes_, system_scalar_fluxes_;

  // Mock pointers and observation pointers
  Framework* mock_framework_obs_ptr_{ nullptr };
  std::shared_ptr<SphericalHarmonicMoments> framework_current_moments_{ nullptr }, system_current_moments_{ nullptr };

  // Test parameters
  const int total_energy_groups_{ test_helpers::RandomInt(5, 10) };
  const int scalar_flux_size_{ test_helpers::RandomInt(10, 20) };

  auto SetUp() -> void override;
};

auto GetScalarFluxFromFrameworkSubroutineTest::SetUp() -> void {
  test_subroutine_ = std::make_unique<GetScalarFluxFromFramework>(std::make_unique<Framework>());
  mock_framework_obs_ptr_ = dynamic_cast<Framework*>(test_subroutine_->framework_ptr());
  framework_system_ptr_ = std::make_shared<System>();

  test_system_.total_groups = total_energy_groups_;

  for (int group = 0; group < total_energy_groups_; ++group) {
    system_scalar_fluxes_.insert({group, Vector(scalar_flux_size_)});
    Vector framework_scalar_flux(scalar_flux_size_);
    auto scalar_flux_values{ test_helpers::RandomVector(scalar_flux_size_, 0, 100) };
    for (int i = 0; i < scalar_flux_size_; ++i) {
      framework_scalar_flux[i] = scalar_flux_values.at(i);
    }
    framework_scalar_fluxes_.insert({group, framework_scalar_flux});
  }

  framework_current_moments_ = std::make_shared<SphericalHarmonicMoments>();
  framework_system_ptr_->current_moments = framework_current_moments_;
  system_current_moments_ = std::make_shared<SphericalHarmonicMoments>();
  test_system_.current_moments = system_current_moments_;
}

TEST_F(GetScalarFluxFromFrameworkSubroutineTest, ConstructorAndGetter) {
  std::unique_ptr<GetScalarFluxFromFramework> test_subroutine{ nullptr };
  EXPECT_NO_THROW({
    test_subroutine = std::make_unique<GetScalarFluxFromFramework>(std::make_unique<Framework>());
  });
  ASSERT_NE(test_subroutine->framework_ptr(), nullptr);
}

TEST_F(GetScalarFluxFromFrameworkSubroutineTest, ConstructorBadDependency) {
  std::unique_ptr<GetScalarFluxFromFramework> test_subroutine{ nullptr };
  const int n_dependencies{ 1 };
  for (int i = 0; i < n_dependencies; ++i) {
    EXPECT_ANY_THROW({
      test_subroutine = std::make_unique<GetScalarFluxFromFramework>(i == 0 ? nullptr : std::make_unique<Framework>());
                     });
  }
}

TEST_F(GetScalarFluxFromFrameworkSubroutineTest, Execute) {
  EXPECT_CALL(*mock_framework_obs_ptr_, SolveSystem());
  EXPECT_CALL(*mock_framework_obs_ptr_, system()).WillOnce(Return(framework_system_ptr_.get()));

  for (int group = 0; group < total_energy_groups_; ++group) {
    std::array<int, 3> index{group, 0, 0};
    EXPECT_CALL(*framework_current_moments_, GetMoment(index)).WillOnce(ReturnRef(framework_scalar_fluxes_.at(group)));
    EXPECT_CALL(*system_current_moments_, GetMoment(index)).WillOnce(ReturnRef(system_scalar_fluxes_.at(group)));
  }

  test_subroutine_->Execute(test_system_);

  for (auto&[group, framework_scalar_flux] : framework_scalar_fluxes_) {
    EXPECT_EQ(framework_scalar_flux, system_scalar_fluxes_.at(group));
  }
}

} // namespace