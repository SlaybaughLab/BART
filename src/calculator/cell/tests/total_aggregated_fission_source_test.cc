#include "calculator/cell/total_aggregated_fission_source.h"

#include <memory>

#include "domain/tests/definition_mock.h"
#include "system/moments/tests/spherical_harmonic_mock.h"
#include "calculator/cell/tests/integrated_fission_source_mock.h"
#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

using ::testing::Return, ::testing::_;

template <typename DimensionWrapper>
class TotalAggregatedFissionSourceTest :
    public ::testing::Test,
    public bart::testing::DealiiTestDomain<DimensionWrapper::value> {
 protected:
  static constexpr int dim = DimensionWrapper::value;

  using CellValueType = calculator::cell::IntegratedFissionSourceMock<dim>;

  // Supporting mock objects
  std::unique_ptr<CellValueType> cell_value_ptr_;
  std::shared_ptr<domain::DefinitionMock<dim>> domain_ptr_;
  std::shared_ptr<system::moments::SphericalHarmonicMock> moments_ptr_;

  void SetUp() override;
};

template <typename DimensionWrapper>
void TotalAggregatedFissionSourceTest<DimensionWrapper>::SetUp() {
  this->SetUpDealii();
  cell_value_ptr_ = std::make_unique<CellValueType>();
  moments_ptr_ = std::make_shared<system::moments::SphericalHarmonicMock>();
  domain_ptr_ = std::make_shared<domain::DefinitionMock<dim>>();
}


TYPED_TEST_CASE(TotalAggregatedFissionSourceTest, bart::testing::AllDimensions);

TYPED_TEST(TotalAggregatedFissionSourceTest, Constructor) {
  //auto& moments_ptr_ = this->moments_ptr_;
  //auto& domain_ptr_ = this->domain_ptr_;
  constexpr int dim = this->dim;

  using CellValueType = calculator::cell::IntegratedFissionSourceI<dim>;

  calculator::cell::TotalAggregatedFissionSource<dim>
      test_aggregator(std::move(this->cell_value_ptr_),
                      this->domain_ptr_);
  EXPECT_EQ(this->domain_ptr_.use_count(), 2);
}

TYPED_TEST(TotalAggregatedFissionSourceTest, AggregatedFissionSourceMPI) {
  auto& moments_ptr_ = this->moments_ptr_;
  auto& domain_ptr_ = this->domain_ptr_;
  constexpr int dim = this->dim;

  calculator::cell::TotalAggregatedFissionSource<dim>
      test_aggregator(std::move(this->cell_value_ptr_),
                      this->domain_ptr_);

  int total_active_cells = this->triangulation_.n_global_active_cells();

  if (dim == 1)
    total_active_cells = dealii::Utilities::MPI::sum(total_active_cells, MPI_COMM_WORLD);

  double cell_value = 5.78;
  double expected_source = total_active_cells * cell_value;

  auto cell_value_obs_ptr =
      dynamic_cast<calculator::cell::IntegratedFissionSourceMock<dim>*>(test_aggregator.cell_fission_source_ptr());

  EXPECT_CALL(*domain_ptr_, Cells())
      .WillOnce(Return(this->cells_));

  for (auto& cell : this->cells_) {
    EXPECT_CALL(*cell_value_obs_ptr, CellValue(cell, moments_ptr_.get()))
        .WillOnce(Return(cell_value));
  }

  double source = test_aggregator.AggreatedFissionSource(moments_ptr_.get());

  EXPECT_NEAR(source, expected_source, 1e-12);


}

} // namespace