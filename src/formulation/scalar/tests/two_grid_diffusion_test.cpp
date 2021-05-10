#include "formulation/scalar/two_grid_diffusion.hpp"

#include "data/cross_sections/tests/cross_sections_mock.hpp"
#include "data/cross_sections/tests/one_group_cross_sections_mock.hpp"
#include "domain/domain_types.hpp"
#include "domain/finite_element/tests/finite_element_mock.hpp"
#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"
#include "test_helpers/test_assertions.hpp"

namespace {

using namespace bart;

using ::testing::NiceMock, ::testing::DoDefault, ::testing::Return, ::testing::_, ::testing::AtLeast;

template <typename DimensionWrapper>
class TwoGridDiffusionTest : public ::testing::Test, public bart::testing::DealiiTestDomain<DimensionWrapper::value> {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using CellPtr = domain::CellPtr<dim>;
  using CrossSectionsMock = data::cross_sections::CrossSectionsMock;
  using OneGroupCrossSectionsMock = data::cross_sections::OneGroupCrossSectionsMock;
  using FiniteElementMock = domain::finite_element::FiniteElementMock<dim>;
  auto SetUp() -> void override;

  // Dependencies
  std::shared_ptr<FiniteElementMock> finite_element_mock_ptr_{ std::make_shared<FiniteElementMock>() };
  std::shared_ptr<CrossSectionsMock> cross_sections_mock_ptr_{ std::make_shared<CrossSectionsMock>() };
  std::shared_ptr<OneGroupCrossSectionsMock> one_group_cross_sections_mock_ptr_{
      std::make_shared<OneGroupCrossSectionsMock>() };

  CellPtr test_cell_ptr_;
  // Test parameters
  static constexpr int cell_dofs{ 2 };
  static constexpr int cell_quad_pts{ 4 };

  const int test_material_id_{ test_helpers::RandomInt(0, 10) };
};

template <typename DimensionWrapper>
auto TwoGridDiffusionTest<DimensionWrapper>::SetUp() -> void {
  this->SetUpDealii();
  test_cell_ptr_ = *this->cells_.begin();
  test_cell_ptr_->set_material_id(test_material_id_);

  ON_CALL(*finite_element_mock_ptr_, dofs_per_cell()).WillByDefault(Return(cell_dofs));
  ON_CALL(*finite_element_mock_ptr_, n_cell_quad_pts()).WillByDefault(Return(cell_quad_pts));
  ON_CALL(*finite_element_mock_ptr_, n_face_quad_pts()).WillByDefault(Return(cell_quad_pts));

  for (int q = 0; q < cell_quad_pts; ++q) {
    ON_CALL(*finite_element_mock_ptr_, Jacobian(q)).WillByDefault(Return((1 + q) * 3));
    for (int i = 0; i < cell_dofs; ++i) {
      ON_CALL(*finite_element_mock_ptr_, ShapeValue(i, q)).WillByDefault(Return(i + q));
      ON_CALL(*finite_element_mock_ptr_, ShapeGradient(i, q)).WillByDefault(Return(dealii::Tensor<1, dim>()));
    }
  }
  ON_CALL(*one_group_cross_sections_mock_ptr_, SigmaAbsorption(test_material_id_)).WillByDefault(Return(0.12));
}

TYPED_TEST_SUITE(TwoGridDiffusionTest, bart::testing::AllDimensions);

TYPED_TEST(TwoGridDiffusionTest, ConstructorAndGettersTest) {
  constexpr int dim{ this->dim };
  EXPECT_CALL(*this->finite_element_mock_ptr_, dofs_per_cell()).WillOnce(DoDefault());
  EXPECT_CALL(*this->finite_element_mock_ptr_, n_cell_quad_pts()).WillOnce(DoDefault());
  EXPECT_CALL(*this->finite_element_mock_ptr_, n_face_quad_pts()).WillOnce(DoDefault());

  formulation::scalar::TwoGridDiffusion<dim> test_formulation(this->finite_element_mock_ptr_,
                                                              this->cross_sections_mock_ptr_,
                                                              this->one_group_cross_sections_mock_ptr_);
  EXPECT_EQ(test_formulation.finite_element_ptr(), this->finite_element_mock_ptr_.get());
  EXPECT_EQ(test_formulation.cross_sections_ptr(), this->cross_sections_mock_ptr_.get());
  EXPECT_EQ(test_formulation.one_group_cross_sections_ptr(), this->one_group_cross_sections_mock_ptr_.get());
}

TYPED_TEST(TwoGridDiffusionTest, FillCellCollisionTerm) {
  using FullMatrix = dealii::FullMatrix<double>;
  constexpr int cell_dofs{ this->cell_dofs };
  constexpr int dim{ this->dim };
  FullMatrix cell_matrix(cell_dofs, cell_dofs);
  const FullMatrix expected_matrix(cell_dofs, cell_dofs, std::vector<double>{18, 25.2, 25.2, 36}.data());

  EXPECT_CALL(*this->finite_element_mock_ptr_, dofs_per_cell()).WillOnce(DoDefault());
  EXPECT_CALL(*this->finite_element_mock_ptr_, n_cell_quad_pts()).WillOnce(DoDefault());
  EXPECT_CALL(*this->finite_element_mock_ptr_, n_face_quad_pts()).WillOnce(DoDefault());

  formulation::scalar::TwoGridDiffusion<dim> test_formulation(this->finite_element_mock_ptr_,
                                                              this->cross_sections_mock_ptr_,
                                                              this->one_group_cross_sections_mock_ptr_);

  EXPECT_CALL(*this->finite_element_mock_ptr_, SetCell(this->test_cell_ptr_)).Times(2); // Once for each call
  EXPECT_CALL(*this->finite_element_mock_ptr_, ShapeValue(_,_))
      .Times(AtLeast(this->cell_dofs * this->cell_quad_pts)) // At least once for each combination of DOF and Quad pt
      .WillRepeatedly(DoDefault());
  EXPECT_CALL(*this->finite_element_mock_ptr_, ShapeGradient(_,_))
      .Times(AtLeast(this->cell_dofs * this->cell_quad_pts))
      .WillRepeatedly(DoDefault());
  EXPECT_CALL(*this->one_group_cross_sections_mock_ptr_, SigmaAbsorption(this->test_material_id_)).WillOnce(DoDefault());
  EXPECT_CALL(*this->finite_element_mock_ptr_, Jacobian(_))
      .Times(AtLeast(this->cell_quad_pts))
      .WillRepeatedly(DoDefault());

  test_formulation.Precalculate(this->test_cell_ptr_);
  test_formulation.FillCellCollisionTerm(cell_matrix, this->test_cell_ptr_, 0);
  EXPECT_TRUE(test_helpers::AreEqual(cell_matrix, expected_matrix));
}

} // namespace
