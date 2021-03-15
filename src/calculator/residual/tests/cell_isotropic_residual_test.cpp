#include "calculator/residual/cell_isotropic_residual.hpp"

#include "data/cross_sections/tests/cross_sections_mock.hpp"
#include "domain/finite_element/tests/finite_element_mock.hpp"
#include "system/moments/tests/spherical_harmonic_mock.h"
#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"
#include "test_helpers/test_assertions.hpp"

namespace  {

using namespace bart;
using ::testing::NiceMock, ::testing::Return, ::testing::DoDefault, ::testing::AtLeast, ::testing::_, ::testing::ReturnRef;

template <typename DimensionWrapper>
class CellIsotropicResidualTest : public bart::testing::DealiiTestDomain<DimensionWrapper::value>, public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using CellIsotropicResidualCalculator = typename calculator::residual::CellIsotropicResidual<dim>;
  using CrossSectionsMock = data::cross_sections::CrossSectionsMock;
  using FiniteElementMock = NiceMock<domain::finite_element::FiniteElementMock<dim>>;
  using FluxMomentsMock = NiceMock<system::moments::SphericalHarmonicMock>;
  using Vector = dealii::Vector<double>;

  std::unique_ptr<CellIsotropicResidualCalculator> test_calculator_;

  // Dependencies
  std::shared_ptr<CrossSectionsMock> cross_sections_mock_ptr_{ std::make_shared<CrossSectionsMock>() };
  std::shared_ptr<FiniteElementMock> finite_element_mock_ptr_{ std::make_shared<FiniteElementMock>() };
  std::shared_ptr<FluxMomentsMock> current_flux_moments_ptr_{ std::make_shared<FluxMomentsMock>() };
  std::shared_ptr<FluxMomentsMock> previous_flux_moments_ptr_{ std::make_shared<FluxMomentsMock>() };

  // Support objects
  std::unordered_map<int, Vector> current_flux_moments_, previous_flux_moments_;

  static constexpr int n_cell_quad_{ 3 };
  static constexpr int n_cell_dofs_{ 2 };
  static constexpr int n_groups_{ 4 };
  static constexpr int global_dofs_{ 10 };
  const int material_id_{ test_helpers::RandomInt(0, 10) };
  auto SetUp() -> void;
};

template <typename DimensionWrapper>
auto CellIsotropicResidualTest<DimensionWrapper>::SetUp() -> void {
  this->SetUpDealii();

  ON_CALL(*current_flux_moments_ptr_, total_groups()).WillByDefault(Return(n_groups_));
  ON_CALL(*previous_flux_moments_ptr_, total_groups()).WillByDefault(Return(n_groups_));

  for (int group = 0; group < n_groups_; ++group) {
    const std::array<int, 3> index{ group, 0, 0};
    auto current_flux_values{ test_helpers::RandomVector(global_dofs_, 0, 100) };
    auto previous_flux_values{ test_helpers::RandomVector(global_dofs_, 0, 100) };

    current_flux_moments_[group] = Vector(current_flux_values.cbegin(), current_flux_values.cend());
    previous_flux_moments_[group] = Vector(previous_flux_values.cbegin(), previous_flux_values.cend());

    ON_CALL(*current_flux_moments_ptr_, GetMoment(index)).WillByDefault(ReturnRef(current_flux_moments_.at(group)));
    ON_CALL(*previous_flux_moments_ptr_, GetMoment(index)).WillByDefault(ReturnRef(previous_flux_moments_.at(group)));
  }

  ON_CALL(*finite_element_mock_ptr_, n_cell_quad_pts()).WillByDefault(Return(n_cell_quad_));
  for (int q = 0; q < n_cell_quad_; ++q) {
    ON_CALL(*finite_element_mock_ptr_, Jacobian(q)).WillByDefault(Return((q + 1) * 3));
  }
  test_calculator_ = std::make_unique<CellIsotropicResidualCalculator>(cross_sections_mock_ptr_,
                                                                       finite_element_mock_ptr_);

  for (auto cell : this->cells_)
    cell->set_material_id(this->material_id_);
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
  using FiniteElementMock = NiceMock<domain::finite_element::FiniteElementMock<this->dim>>;
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
  const auto& test_cell = *this->cells_.begin();
  constexpr int test_group{ 1 };
  EXPECT_CALL(*this->finite_element_mock_ptr_, SetCell(test_cell));
  std::unordered_map<int, dealii::FullMatrix<double>> sigma_s;
  sigma_s[this->material_id_] = dealii::FullMatrix<double>(
      this->n_groups_, this->n_groups_, std::array<double, 16>{0.5, 0, 0, 0,
                                                               0.2, 0.7, 0.01, 0.1,
                                                               0.3, 0.5, 0.6, 0.1,
                                                               0.4, 0.6, 0.9, 0.8}.data());
  EXPECT_CALL(*this->cross_sections_mock_ptr_, sigma_s()).WillOnce(Return(sigma_s));

  EXPECT_CALL(*this->finite_element_mock_ptr_, ValueAtQuadrature(this->current_flux_moments_.at(2)))
      .WillOnce(Return(std::vector<double>{3, 3, 3}));
  EXPECT_CALL(*this->finite_element_mock_ptr_, ValueAtQuadrature(this->previous_flux_moments_.at(2)))
      .WillOnce(Return(std::vector<double>{1, 1, 1}));
  EXPECT_CALL(*this->finite_element_mock_ptr_, ValueAtQuadrature(this->current_flux_moments_.at(3)))
      .WillOnce(Return(std::vector<double>{5, 5, 5}));
  EXPECT_CALL(*this->finite_element_mock_ptr_, ValueAtQuadrature(this->previous_flux_moments_.at(3)))
      .WillOnce(Return(std::vector<double>{2, 2, 2}));

  const double expected_result{ 5.76 };

  const auto result = this->test_calculator_->CalculateCellResidual(test_cell, this->current_flux_moments_ptr_.get(),
                                                                    this->previous_flux_moments_ptr_.get(),
                                                                    test_group);
  EXPECT_NEAR(expected_result, result, 1e-10);
}

} // namespace

