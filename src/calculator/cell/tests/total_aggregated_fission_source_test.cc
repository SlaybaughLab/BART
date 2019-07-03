#include "calculator/cell/total_aggregated_fission_source.h"

#include <memory>

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

  std::unique_ptr<CellValueType> cell_value_ptr_;
  std::shared_ptr<system::moments::SphericalHarmonicMock> moments_ptr_;

  void SetUp() override;
};

template <typename DimensionWrapper>
void TotalAggregatedFissionSourceTest<DimensionWrapper>::SetUp() {
  this->SetUpDealii();
  cell_value_ptr_ = std::make_unique<CellValueType>();
  moments_ptr_ = std::make_shared<system::moments::SphericalHarmonicMock>();
}


TYPED_TEST_CASE(TotalAggregatedFissionSourceTest, bart::testing::AllDimensions);

TYPED_TEST(TotalAggregatedFissionSourceTest, Constructor) {
  auto& moments_ptr_ = this->moments_ptr_;
  constexpr int dim = this->dim;

  using CellValueType = calculator::cell::IntegratedFissionSourceI<dim>;

  calculator::cell::TotalAggregatedFissionSource<dim>
      test_aggregator(std::move(this->cell_value_ptr_));

}

} // namespace