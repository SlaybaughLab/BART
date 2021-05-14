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

  static constexpr int n_groups_{ 3 };
  const int material_id_{ test_helpers::RandomInt(0, 10) };
  auto SetUp() -> void;
};

template <typename DimensionWrapper>
auto CellIsotropicResidualTest<DimensionWrapper>::SetUp() -> void {
  this->SetUpDealii();

  ON_CALL(*current_flux_moments_ptr_, total_groups()).WillByDefault(Return(n_groups_));
  ON_CALL(*previous_flux_moments_ptr_, total_groups()).WillByDefault(Return(n_groups_));
  ON_CALL(*finite_element_mock_ptr_, dofs_per_cell()).WillByDefault(Return(this->fe_.dofs_per_cell));

  for (int group = 0; group < n_groups_; ++group) {
    const std::array<int, 3> index{ group, 0, 0};

    current_flux_moments_[group] = Vector(this->dof_handler_.n_dofs());
    previous_flux_moments_[group] = Vector(this->dof_handler_.n_dofs());
    current_flux_moments_[group].add(10 * (group + 1));
    previous_flux_moments_[group].add(group + 1);

    ON_CALL(*current_flux_moments_ptr_, GetMoment(index)).WillByDefault(ReturnRef(current_flux_moments_.at(group)));
    ON_CALL(*previous_flux_moments_ptr_, GetMoment(index)).WillByDefault(ReturnRef(previous_flux_moments_.at(group)));
  }

  using MaterialIDMappedToSigmaS = CrossSectionsMock::MaterialIDMappedTo<dealii::FullMatrix<double>>;
  std::array<double, 9> sigma_s_values{2, 0.5, 0.25, 1.0/3.0, 3, 2.0/3.0, 1.0/5.0, 2.0/5.0, 5};
  dealii::FullMatrix<double> sigma_s_matrix(3, 3, sigma_s_values.data());
  MaterialIDMappedToSigmaS sigma_s_map{{material_id_, sigma_s_matrix}};
  ON_CALL(*cross_sections_mock_ptr_, sigma_s()).WillByDefault(Return(sigma_s_map));

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

TYPED_TEST(CellIsotropicResidualTest, CalculateCellResidual) {
  const int n_cells = this->cells_.size();
  const int total_calls{ n_cells * this->n_groups_ };

  EXPECT_CALL(*this->current_flux_moments_ptr_, total_groups()).Times(total_calls).WillRepeatedly(DoDefault());
  EXPECT_CALL(*this->cross_sections_mock_ptr_, sigma_s()).Times(total_calls).WillRepeatedly(DoDefault());
  EXPECT_CALL(*this->finite_element_mock_ptr_, dofs_per_cell()).Times(total_calls).WillRepeatedly(DoDefault());



  dealii::Vector<double> expected_residual(this->dof_handler_.n_dofs());
  std::fill(expected_residual.begin(), expected_residual.end(), 33.75);
  dealii::Vector<double> calculated_residual(this->dof_handler_.n_dofs());

  for (int group = 0; group < this->n_groups_; ++group) {
    for (int group_in = group + 1; group_in < this->n_groups_; ++group_in) {
      std::array<int, 3> index{group_in, 0, 0};
      EXPECT_CALL(*this->current_flux_moments_ptr_, GetMoment(index))
          .Times(AtLeast(1)).WillRepeatedly(DoDefault());
      EXPECT_CALL(*this->previous_flux_moments_ptr_, GetMoment(index))
          .Times(AtLeast(1)).WillRepeatedly(DoDefault());
    }
    for (const auto cell : this->cells_) {
      this->test_calculator_->CalculateCellResidual(calculated_residual, cell,
                                                    this->current_flux_moments_ptr_.get(),
                                                    this->previous_flux_moments_ptr_.get(),
                                                    group);
    }
  }
  for (int i = 0; i < this->dof_handler_.n_dofs(); ++i) {
    EXPECT_NEAR(std::fmod(calculated_residual(i), 33.75), 0, 1e-10);
  }
}

} // namespace

