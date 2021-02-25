#include "acceleration/two_grid/flux_corrector.hpp"

#include "acceleration/two_grid/tests/flux_corrector_mock.hpp"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

using namespace bart;
using ::testing::ContainerEq;

class AccelerationTwoGridFluxCorrectorTest : public ::testing::Test {
 public:
  using GroupToDomainSpectralShapeMap = std::unordered_map<int, dealii::Vector<double>>;
  using TestFluxCorrector = acceleration::two_grid::FluxCorrector;
  using Vector = dealii::Vector<double>;

  // Test supporting objects
  Vector flux, error, corrected_flux;
  GroupToDomainSpectralShapeMap group_to_domain_spectral_shape_map_;

  // Test parameters
  const int total_groups{ test_helpers::RandomInt(5, 10) };
  const int total_dofs{ test_helpers::RandomInt(10, 20) };
  const int test_group{ test_helpers::RandomInt(0, total_groups - 1)};

  auto SetUp() -> void override;
};

auto AccelerationTwoGridFluxCorrectorTest::SetUp() -> void {
  for (int group = 0; group < total_groups; ++group) {
    auto spectral_shape_vector = test_helpers::RandomVector(total_dofs, 0, 1);
    group_to_domain_spectral_shape_map_[group] = dealii::Vector<double>(spectral_shape_vector.begin(),
                                                                        spectral_shape_vector.end());
  }
  const auto flux_vector = test_helpers::RandomVector(total_dofs, 0, 100);
  const auto error_vector = test_helpers::RandomVector(total_dofs, 0, 100);
  const auto spectral_shape = group_to_domain_spectral_shape_map_.at(test_group);

  flux.reinit(total_dofs);
  error.reinit(total_dofs);
  corrected_flux.reinit(total_dofs);

  for (int i = 0; i < total_dofs; ++i) {
    flux[i] = flux_vector.at(i);
    error[i] = error_vector.at(i);
    corrected_flux[i] = flux_vector.at(i) + spectral_shape[i] * error_vector.at(i);
  }
}

TEST_F(AccelerationTwoGridFluxCorrectorTest, ConstructorAndGetter) {
  TestFluxCorrector test_corrector(this->group_to_domain_spectral_shape_map_);
  EXPECT_THAT(test_corrector.group_to_domain_spectral_shape_map(), ContainerEq(group_to_domain_spectral_shape_map_));
}

TEST_F(AccelerationTwoGridFluxCorrectorTest, CorrectFlux) {
  TestFluxCorrector test_corrector(this->group_to_domain_spectral_shape_map_);
  test_corrector.CorrectFlux(flux, error, test_group);
  EXPECT_EQ(flux, corrected_flux);
}


} // namespace
