#include "calculator/cell/total_aggregated_fission_source.hpp"

#include <memory>

#include "domain/tests/domain_mock.h"
#include "system/moments/tests/spherical_harmonic_mock.h"
#include "calculator/cell/tests/integrated_fission_source_mock.hpp"
#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

using ::testing::Return, ::testing::_;

template <typename DimensionWrapper>
class TotalAggregatedFissionSourceTest : public ::testing::Test,
                                         public bart::testing::DealiiTestDomain<DimensionWrapper::value> {
 public:
  static constexpr int dim = DimensionWrapper::value;

  using TestAggregatedFissionSource = calculator::cell::TotalAggregatedFissionSource<dim>;
  using DomainMock = domain::DefinitionMock<dim>;
  using IntegratedFissionSourceMock = calculator::cell::IntegratedFissionSourceMock<dim>;
  using MomentsMock = system::moments::SphericalHarmonicMock;

  // Constructed test object
  std::unique_ptr<TestAggregatedFissionSource> test_aggregated_fission_source_;

  // Supporting mock objects
  std::shared_ptr<DomainMock> domain_ptr_{ std::make_shared<DomainMock>() };
  std::shared_ptr<MomentsMock> moments_ptr_{ std::make_shared<MomentsMock>() };

  IntegratedFissionSourceMock* cell_integrated_fission_source_obs_ptr_{ nullptr };

  void SetUp() override;
};

template <typename DimensionWrapper>
void TotalAggregatedFissionSourceTest<DimensionWrapper>::SetUp() {
  this->SetUpDealii();
  auto integrated_fission_source_mock_ptr = std::make_unique<IntegratedFissionSourceMock>();
  cell_integrated_fission_source_obs_ptr_ = integrated_fission_source_mock_ptr.get();

  test_aggregated_fission_source_ = std::make_unique<TestAggregatedFissionSource>(
      std::move(integrated_fission_source_mock_ptr), domain_ptr_);
}

TYPED_TEST_CASE(TotalAggregatedFissionSourceTest, bart::testing::AllDimensions);

/* Constructor should not throw on good dependencies */
TYPED_TEST(TotalAggregatedFissionSourceTest, Constructor) {
  constexpr int dim = this->dim;
  using TestAggregatedFissionSource = calculator::cell::TotalAggregatedFissionSource<dim>;
  using DomainMock = domain::DefinitionMock<dim>;
  using IntegratedFissionSourceMock = calculator::cell::IntegratedFissionSourceMock<dim>;

  EXPECT_NO_THROW({
    TestAggregatedFissionSource test_fission_source(std::make_unique<IntegratedFissionSourceMock>(),
                                                    std::make_shared<DomainMock>());
  });
}

/* Constructor should throw on any null dependencies */
TYPED_TEST(TotalAggregatedFissionSourceTest, ConstructorBadDependencies) {
  constexpr int dim = this->dim;
  using TestAggregatedFissionSource = calculator::cell::TotalAggregatedFissionSource<dim>;
  using DomainMock = domain::DefinitionMock<dim>;
  using IntegratedFissionSourceMock = calculator::cell::IntegratedFissionSourceMock<dim>;

  EXPECT_ANY_THROW({
    TestAggregatedFissionSource test_fission_source(nullptr, std::make_shared<DomainMock>());
  });
  EXPECT_ANY_THROW({
    TestAggregatedFissionSource test_fission_source(std::make_unique<IntegratedFissionSourceMock>(), nullptr);
  });
}

/* Getters should return the correct pointers to the dependencies */
TYPED_TEST(TotalAggregatedFissionSourceTest, Getters) {
  EXPECT_EQ(this->test_aggregated_fission_source_->cell_fission_source_ptr(),
            this->cell_integrated_fission_source_obs_ptr_);
  EXPECT_EQ(this->test_aggregated_fission_source_->domain_ptr(), this->domain_ptr_.get());
}

/* Calculated aggregated fission source should be properly calculated. This is done by mocking the return value of the
 * cell integrator and ensuring this class is adding it all properly AND that this works with MPI. */
TYPED_TEST(TotalAggregatedFissionSourceTest, AggregatedFissionSourceMPI) {
  auto& moments_ptr_ = this->moments_ptr_;
  auto& domain_ptr_ = this->domain_ptr_;
  auto& test_aggregator = *this->test_aggregated_fission_source_;
  constexpr int dim = this->dim;

  // Get total active cells, an extra step is required if 1D due to the way triangulations work in Dealii.
  int total_active_cells = this->triangulation_.n_global_active_cells();
  if (dim == 1)
    total_active_cells = dealii::Utilities::MPI::sum(total_active_cells, MPI_COMM_WORLD);

  const double cell_value{ 5.78 }, expected_source{ total_active_cells * cell_value };

  EXPECT_CALL(*domain_ptr_, Cells()).WillOnce(Return(this->cells_));
  for (auto& cell : this->cells_) {
    EXPECT_CALL(*this->cell_integrated_fission_source_obs_ptr_, CellValue(cell, moments_ptr_.get()))
        .WillOnce(Return(cell_value));
  }

  const double source = test_aggregator.AggregatedFissionSource(moments_ptr_.get());
  EXPECT_NEAR(source, expected_source, 1e-12);
}

} // namespace