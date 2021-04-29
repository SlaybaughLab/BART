#include "calculator/residual/tests/cell_isotropic_residual_mock.hpp"
#include "calculator/residual/domain_isotropic_residual.hpp"
#include "domain/tests/domain_mock.hpp"
#include "system/moments/tests/spherical_harmonic_mock.h"
#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_assertions.hpp"

namespace  {

using namespace bart;

using ::testing::DoDefault, ::testing::NiceMock, ::testing::Return, ::testing::_, ::testing::AtLeast;

template <typename DimensionWrapper>
class CalculatorResidualDomainIsotropicResidualTest : public ::testing::Test,
                                                      public bart::testing::DealiiTestDomain<DimensionWrapper::value> {
 public:
  static constexpr int dim{ DimensionWrapper::value };
  using CellIsotropicResidualMock = typename calculator::residual::CellIsotropicResidualMock<dim>;
  using DomainMock = typename domain::DomainMock<dim>;
  using FluxMomentsMock = system::moments::SphericalHarmonicMock;
  using Vector = dealii::Vector<double>;

  using TestDomainIsotropicResidual = typename calculator::residual::DomainIsotropicResidual<dim>;

  // Test object
  std::unique_ptr<TestDomainIsotropicResidual> test_calculator_;

  // Dependencies
  CellIsotropicResidualMock* cell_isotropic_residual_mock_obs_ptr_;
  std::shared_ptr<DomainMock> domain_mock_ptr_{ std::make_shared<DomainMock>() };

  // Test parameters and supporting objects
  Vector expected_isotropic_residual_;
  std::shared_ptr<FluxMomentsMock> current_flux_moments_{ std::make_shared<FluxMomentsMock>() };
  std::shared_ptr<FluxMomentsMock> previous_flux_moments_{ std::make_shared<FluxMomentsMock>() };

  static constexpr int total_groups{ 3 };
  auto SetUp() -> void override;
};

template <typename DimensionWrapper>
auto CalculatorResidualDomainIsotropicResidualTest<DimensionWrapper>::SetUp() -> void {
  this->SetUpDealii();

  ON_CALL(*current_flux_moments_, total_groups()).WillByDefault(Return(total_groups));
  ON_CALL(*previous_flux_moments_, total_groups()).WillByDefault(Return(total_groups));

  auto cell_isotropic_residual_mock_ptr = std::make_unique<CellIsotropicResidualMock>();
  cell_isotropic_residual_mock_obs_ptr_ = cell_isotropic_residual_mock_ptr.get();

  test_calculator_ = std::make_unique<TestDomainIsotropicResidual>(std::move(cell_isotropic_residual_mock_ptr),
                                                                   domain_mock_ptr_);
}

TYPED_TEST_SUITE(CalculatorResidualDomainIsotropicResidualTest, bart::testing::AllDimensions);

TYPED_TEST(CalculatorResidualDomainIsotropicResidualTest, DependencyGetters) {
  EXPECT_EQ(this->test_calculator_->domain_ptr(), this->domain_mock_ptr_.get());
  EXPECT_EQ(this->test_calculator_->cell_isotropic_residual_calculator_ptr(),
            this->cell_isotropic_residual_mock_obs_ptr_);
}

TYPED_TEST(CalculatorResidualDomainIsotropicResidualTest, NullDependenciesThrow) {
  constexpr int n_dependencies{ 2 };
  constexpr int dim{ this->dim };
  using CellIsotropicResidualMock = typename calculator::residual::CellIsotropicResidualMock<dim>;
  using DomainMock = typename domain::DomainMock<dim>;
  using TestDomainIsotropicResidual = typename calculator::residual::DomainIsotropicResidual<dim>;
  for (int i = 0; i < n_dependencies; ++i) {
    EXPECT_ANY_THROW({
      TestDomainIsotropicResidual(i == 0 ? nullptr : std::make_unique<CellIsotropicResidualMock>(),
                                  i == 1 ? nullptr : std::make_shared<DomainMock>());
    });
  }
}

TYPED_TEST(CalculatorResidualDomainIsotropicResidualTest, CalculateDomainResidualBadTotalGroups) {
  EXPECT_CALL(*this->current_flux_moments_, total_groups()).WillOnce(DoDefault());
  EXPECT_CALL(*this->previous_flux_moments_, total_groups()).WillOnce(Return(this->total_groups + 1));
  EXPECT_ANY_THROW({
    this->test_calculator_->CalculateDomainResidual(this->current_flux_moments_.get(),
                                                    this->previous_flux_moments_.get());
  });
}

TYPED_TEST(CalculatorResidualDomainIsotropicResidualTest, CalculateDomainResidual) {
  EXPECT_CALL(*this->current_flux_moments_, total_groups()).WillOnce(DoDefault());
  EXPECT_CALL(*this->previous_flux_moments_, total_groups()).WillOnce(DoDefault());
  EXPECT_CALL(*this->domain_mock_ptr_, Cells()).Times(AtLeast(1)).WillRepeatedly(Return(this->cells_));
  EXPECT_CALL(*this->domain_mock_ptr_, total_degrees_of_freedom())
      .Times(AtLeast(1))
      .WillRepeatedly(Return(this->dof_handler_.n_dofs()));
  EXPECT_CALL(*this->domain_mock_ptr_, GetCellVector())
      .WillOnce(Return(dealii::Vector<double>(this->fe_.dofs_per_cell)));

  for (int group = 0; group < this->total_groups; ++group) {
    for (auto& cell : this->cells_) {
      EXPECT_CALL(*this->cell_isotropic_residual_mock_obs_ptr_,
                  CalculateCellResidual(_, cell, this->current_flux_moments_.get(),
                                        this->previous_flux_moments_.get(), group));
    }
  }

  const auto result = this->test_calculator_->CalculateDomainResidual(this->current_flux_moments_.get(),
                                                                      this->previous_flux_moments_.get());
  ASSERT_EQ(result.size(), this->dof_handler_.n_dofs());
}

} // namespace
