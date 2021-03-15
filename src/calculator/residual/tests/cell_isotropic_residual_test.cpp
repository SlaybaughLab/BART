#include "calculator/residual/cell_isotropic_residual.hpp"

#include "data/cross_sections/tests/cross_sections_mock.hpp"
#include "domain/finite_element/tests/finite_element_mock.hpp"
#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;
using ::testing::NiceMock, ::testing::Return, ::testing::DoDefault, ::testing::AtLeast, ::testing::_;

template <typename DimensionWrapper>
class CellIsotropicResidualTest : public bart::testing::DealiiTestDomain<DimensionWrapper::value>, public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using CellIsotropicResidualCalculator = typename calculator::residual::CellIsotropicResidual<dim>;
  using CrossSectionsMock = data::cross_sections::CrossSectionsMock;
  using FiniteElementMock = domain::finite_element::FiniteElementMock<dim>;

  std::unique_ptr<CellIsotropicResidualCalculator> test_calculator_;

  // Dependencies
  std::shared_ptr<CrossSectionsMock> cross_sections_mock_ptr_{ std::make_shared<CrossSectionsMock>() };
  std::shared_ptr<FiniteElementMock> finite_element_mock_ptr_{ std::make_shared<FiniteElementMock>() };



  static constexpr int n_cell_quad_{ 3 };
  static constexpr int n_cell_dofs_{ 2 };
  auto SetUp() -> void;
};

template <typename DimensionWrapper>
auto CellIsotropicResidualTest<DimensionWrapper>::SetUp() -> void {
  this->SetUpDealii();
  test_calculator_ = std::make_unique<CellIsotropicResidualCalculator>(cross_sections_mock_ptr_,
                                                                       finite_element_mock_ptr_);

  for (int q = 0; q < n_cell_dofs_; ++q) {
    ON_CALL(*finite_element_mock_ptr_, Jacobian(q)).WillByDefault(Return((q + 1) * 3));
  }
}

TYPED_TEST_SUITE(CellIsotropicResidualTest, bart::testing::AllDimensions);

// Getters should return correct dependency pointers.
TYPED_TEST(CellIsotropicResidualTest, DependencyGetters) {
  EXPECT_EQ(this->test_calculator_->cross_sections_ptr(), this->cross_sections_mock_ptr_.get());
  EXPECT_EQ(this->test_calculator_->finite_element_ptr(), this->finite_element_mock_ptr_.get());
}

// Null dependencies given to the constructor should throw.
TYPED_TEST(CellIsotropicResidualTest, NullDependenciesThrow) {
  using CellIsotropicResidualCalculator = typename calculator::residual::CellIsotropicResidual<this->dim>;
  using CrossSectionsMock = data::cross_sections::CrossSectionsMock;
  using FiniteElementMock = domain::finite_element::FiniteElementMock<this->dim>;
  constexpr int n_dependencies { 2 };
  for (int i = 0; i < n_dependencies; ++i) {
    EXPECT_ANY_THROW({
                       CellIsotropicResidualCalculator( i == 0 ? nullptr : std::make_shared<CrossSectionsMock>(),
                                                        i == 1 ? nullptr : std::make_shared<FiniteElementMock>());
                     });
  }
}

// Call to CalculateCellResidual should calculate correct value
TYPED_TEST(CellIsotropicResidualTest, CalculateCellResidual) {

}

} // namespace

