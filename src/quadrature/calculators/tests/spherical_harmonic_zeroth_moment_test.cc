#include "quadrature/calculators/spherical_harmonic_zeroth_moment.h"

#include <memory>

#include "quadrature/angular/tests/angular_quadrature_set_mock.h"
#include "system/solution/tests/mpi_angular_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

using ::testing::Ref;

template <typename DimensionWrapper>
class QuadCalcSphericalHarmonicMomentsOnlyScalar : public ::testing::Test {
 protected:
  static constexpr int dim = DimensionWrapper::value;
  // Aliases
  using AngularQuadratureSetType = quadrature::angular::AngularQuadratureSetMock<dim>;
  using MomentCalculatorType = quadrature::calculators::SphericalHarmonicZerothMoment<dim>;

  // Pointer to tested object
  std::unique_ptr<MomentCalculatorType> test_calculator;

  // Supporting objects
  system::solution::MPIAngularMock mock_solution_;
  std::shared_ptr<AngularQuadratureSetType> mock_angular_quad_;

  // Observing pointers
  AngularQuadratureSetType* angular_quad_obs_ptr_;

  void SetUp() override;
};

template <typename DimensionWrapper>
void QuadCalcSphericalHarmonicMomentsOnlyScalar<DimensionWrapper>::SetUp() {

  // Instantiate mock objects
  mock_angular_quad_ = std::make_shared<AngularQuadratureSetType>();

  // Instantiate object to be tested
  test_calculator = std::make_unique<MomentCalculatorType>(mock_angular_quad_);

  // Set up observation pointers
  angular_quad_obs_ptr_ = dynamic_cast<AngularQuadratureSetType*>(
      test_calculator->angular_quadrature_set_ptr());

}

TYPED_TEST_CASE(QuadCalcSphericalHarmonicMomentsOnlyScalar,
                bart::testing::AllDimensions);

TYPED_TEST(QuadCalcSphericalHarmonicMomentsOnlyScalar, Constructor) {
  auto& mock_angular_quad_ = this->mock_angular_quad_;
  auto& angular_quad_obs_ptr_ = this->angular_quad_obs_ptr_;

  EXPECT_EQ(mock_angular_quad_.use_count(), 2);
  EXPECT_THAT(*mock_angular_quad_, Ref(*angular_quad_obs_ptr_));
}

} // namespace
