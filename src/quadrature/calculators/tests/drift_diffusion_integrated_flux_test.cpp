#include "quadrature/calculators/drift_diffusion_integrated_flux.hpp"

#include <memory>

#include "quadrature/tests/quadrature_set_mock.h"
#include "quadrature/tests/quadrature_point_mock.h"
#include "quadrature/calculators/tests/drift_diffusion_integrated_flux_mock.hpp"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_assertions.hpp"

namespace  {

using namespace bart;

using ::testing::NiceMock, ::testing::Return, ::testing::DoDefault, ::testing::ContainerEq;

template <typename DimensionWrapper>
class DriftDiffusionIntegratedFluxTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using QuadraturePointType = NiceMock<typename quadrature::QuadraturePointMock<dim>>;
  using QuadratureSetType = NiceMock<typename quadrature::QuadratureSetMock<dim>>;
  using TestIntegrator = typename quadrature::calculators::DriftDiffusionIntegratedFlux<dim>;
  using Vector = dealii::Vector<double>;
  using VectorPtr = std::shared_ptr<dealii::Vector<double>>;
  using VectorMap = std::map<quadrature::QuadraturePointIndex, VectorPtr>;

  DriftDiffusionIntegratedFluxTest()
      : expected_result_(expected_result_values_.cbegin(), expected_result_values_.cend()) {}

  // Test parameters
  static constexpr int n_quadrature_points{ 3 };
  static constexpr int n_total_dofs{ 2 };

  // Test object
  std::unique_ptr<TestIntegrator> test_integrator_{ nullptr };

  // Mock objects
  std::shared_ptr<QuadratureSetType> quadrature_set_ptr_{ nullptr };
  std::array<std::shared_ptr<QuadraturePointType>, n_quadrature_points> mock_quadrature_points_;

  // Supporting objects
  VectorMap angular_flux_map_{};
  const std::vector<double> expected_result_values_{5501.5 * dim, 11003 * dim};
  const Vector expected_result_;

  auto SetUp() -> void override;
};

template <typename DimensionWrapper>
auto DriftDiffusionIntegratedFluxTest<DimensionWrapper>::SetUp() -> void {
  quadrature_set_ptr_ = std::make_shared<QuadratureSetType>();
  test_integrator_ = std::move(std::make_unique<TestIntegrator>(quadrature_set_ptr_));

  using Index = quadrature::QuadraturePointIndex;
  const std::array<double, n_quadrature_points> weights{0.15, 0.25, 0.6};
  for (int i = 0; i < n_quadrature_points; ++i) {
    auto quadrature_point_ptr = std::make_shared<QuadraturePointType>();
    mock_quadrature_points_.at(i) = quadrature_point_ptr;
    dealii::Tensor<1, dim> position_tensor;
    for (int j = 0; j < dim; ++j)
      position_tensor[j] = i + 1;
    ON_CALL(*quadrature_point_ptr, weight()).WillByDefault(Return(weights.at(i)));
    ON_CALL(*quadrature_point_ptr, cartesian_position_tensor()).WillByDefault(Return(position_tensor));
    ON_CALL(*quadrature_set_ptr_, GetQuadraturePoint(Index(i))).WillByDefault(Return(quadrature_point_ptr));
    auto vector_ptr = std::make_shared<Vector>(n_total_dofs);
    for (int j = 0; j < n_total_dofs; ++j)
      (*vector_ptr)[j] = std::pow(10.0, i + 1) * (j + 1);
    angular_flux_map_.insert({Index(i), vector_ptr});
  }
  ON_CALL(*quadrature_set_ptr_, size()).WillByDefault(Return(n_quadrature_points));
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

TYPED_TEST(DriftDiffusionIntegratedFluxTest, Integrate) {
  EXPECT_CALL(*this->quadrature_set_ptr_, size()).WillOnce(DoDefault());
  for (int i = 0; i < this->n_quadrature_points; ++i) {
    using Index = quadrature::QuadraturePointIndex;
    EXPECT_CALL(*this->quadrature_set_ptr_, GetQuadraturePoint(Index(i))).WillOnce(DoDefault());
  }
  for (auto& quadrature_point : this->mock_quadrature_points_) {
    EXPECT_CALL(*quadrature_point, weight()).WillOnce(DoDefault());
    EXPECT_CALL(*quadrature_point, cartesian_position_tensor()).WillOnce(DoDefault());
  }
  auto result = this->test_integrator_->Integrate(this->angular_flux_map_);
  EXPECT_TRUE(test_helpers::AreEqual(result, this->expected_result_));
}

TYPED_TEST(DriftDiffusionIntegratedFluxTest, IntegrateBadAngularFluxSize) {
  using Vector = dealii::Vector<double>;
  using VectorPtr = std::shared_ptr<dealii::Vector<double>>;
  using VectorMap = std::map<quadrature::QuadraturePointIndex, VectorPtr>;
  using Index = quadrature::QuadraturePointIndex;

  EXPECT_CALL(*this->quadrature_set_ptr_, size()).WillOnce(DoDefault());

  VectorMap bad_angular_flux_map;
  bad_angular_flux_map.insert({Index(0), std::make_shared<Vector>()});
  EXPECT_ANY_THROW({
    [[maybe_unused]] auto result = this->test_integrator_->Integrate(bad_angular_flux_map);
  });
}





} // namespace
