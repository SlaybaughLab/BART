#include "calculator/residual/tests/cell_isotropic_residual_mock.hpp"
#include "calculator/residual/domain_isotropic_residual.hpp"
#include "domain/tests/domain_mock.hpp"
#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_assertions.hpp"

namespace  {

using namespace bart;

using ::testing::DoDefault, ::testing::NiceMock, ::testing::Return, ::testing::_;

template <typename DimensionWrapper>
class CalculatorResidualDomainIsotropicResidualTest : public ::testing::Test,
                                                      public bart::testing::DealiiTestDomain<DimensionWrapper::value> {
 public:
  static constexpr int dim{ DimensionWrapper::value };
  using CellIsotropicResidualMock = typename calculator::residual::CellIsotropicResidualMock<dim>;
  using DomainMock = typename domain::DomainMock<dim>;
  using Vector = dealii::Vector<double>;

  using TestDomainIsotropicResidual = typename calculator::residual::DomainIsotropicResidual<dim>;

  // Test object
  std::unique_ptr<TestDomainIsotropicResidual> test_calculator_;

  // Dependencies
  CellIsotropicResidualMock* cell_isotropic_residual_mock_obs_ptr_;
  std::shared_ptr<DomainMock> domain_mock_ptr_{ std::make_shared<DomainMock>() };

  // Test parameters and supporting objects
  Vector expected_isotropic_residual_;

  static constexpr int total_groups{ 3 };
  auto SetUp() -> void override;
};

template <typename DimensionWrapper>
auto CalculatorResidualDomainIsotropicResidualTest<DimensionWrapper>::SetUp() -> void {
  this->SetUpDealii();

  expected_isotropic_residual_.reinit(this->dof_handler_.n_dofs());
  for (auto& cell : this->cells_) {
    std::vector<unsigned int> global_dofs(this->fe_.dofs_per_cell);
    cell->get_dof_indices(global_dofs);
    for (int index : global_dofs) {
      expected_isotropic_residual_[index] += total_groups; // Each group should contribute 1
    }
  }

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

} // namespace
