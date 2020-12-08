#include "quadrature/calculators/drift_diffusion_integrated_flux.hpp"

#include <memory>

#include "quadrature/tests/quadrature_set_mock.h"
#include "quadrature/calculators/tests/drift_diffusion_integrated_flux_mock.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class DriftDiffusionIntegratedFluxTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using QuadratureSetType = typename quadrature::QuadratureSetMock<dim>;

  // Mock objects
  std::shared_ptr<QuadratureSetType> quadrature_set_ptr_;

  auto SetUp() -> void override;
};

template <typename DimensionWrapper>
auto DriftDiffusionIntegratedFluxTest<DimensionWrapper>::SetUp() -> void {
  quadrature_set_ptr_ = std::make_shared<QuadratureSetType>();
}

TYPED_TEST_SUITE(DriftDiffusionIntegratedFluxTest, bart::testing::AllDimensions);

TYPED_TEST(DriftDiffusionIntegratedFluxTest, Constructor) {
  constexpr int dim = this->dim;
  using TestIntegrator = typename quadrature::calculators::DriftDiffusionIntegratedFlux<dim>;
  using QuadratureSetType = typename quadrature::QuadratureSetMock<dim>;

  auto quadrature_set_ptr = std::make_shared<QuadratureSetType>();
  TestIntegrator integrator(quadrature_set_ptr);
  ASSERT_NE(nullptr, integrator.quadrature_set_ptr());
  EXPECT_EQ(integrator.quadrature_set_ptr(), quadrature_set_ptr.get());
}

TYPED_TEST(DriftDiffusionIntegratedFluxTest, ConstructorBadDependencies) {
  constexpr int dim = this->dim;
  using TestIntegrator = typename quadrature::calculators::DriftDiffusionIntegratedFlux<dim>;

  EXPECT_ANY_THROW({
                     TestIntegrator integrator(nullptr);
  });
}



} // namespace
