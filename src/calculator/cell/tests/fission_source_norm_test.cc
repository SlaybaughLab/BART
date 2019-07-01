#include "calculator/cell/fission_source_norm.h"

#include <memory>

#include "data/cross_sections.h"
#include "domain/domain_types.h"
#include "domain/tests/finite_element_mock.h"
#include "material/tests/mock_material.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/dealii_test_domain.h"


namespace  {

using namespace bart;

using ::testing::Return, ::testing::DoDefault, ::testing::NiceMock;

template <typename DimensionWrapper>
class CalcCellFissionSourceNormTest :
    public ::testing::Test,
    public bart::testing::DealiiTestDomain<DimensionWrapper::value> {
 protected:
  static constexpr int dim = DimensionWrapper::value;
  using FiniteElementType = typename domain::FiniteElementMock<dim>;
  using FissionSourceNormType = typename calculator::cell::FissionSourceNorm<dim>;

  // Supporting objects and mocks
  std::shared_ptr<NiceMock<FiniteElementType>> finite_element_ptr_;
  std::shared_ptr<data::CrossSections> cross_sections_ptr_;


  NiceMock<btest::MockMaterial> mock_material_;

  void SetUp() override;
};

template <typename DimensionWrapper>
void CalcCellFissionSourceNormTest<DimensionWrapper>::SetUp() {

  std::unordered_map<int, bool> fissile_id_map{{0, true}, {1, false}};

  ON_CALL(mock_material_, GetFissileIDMap())
      .WillByDefault(Return(fissile_id_map));

  finite_element_ptr_ = std::make_shared<NiceMock<FiniteElementType>>();
  cross_sections_ptr_ = std::make_shared<data::CrossSections>(mock_material_);

  ON_CALL(*finite_element_ptr_, n_cell_quad_pts())
      .WillByDefault(Return(4));

  this->SetUpDealii();
}

TYPED_TEST_CASE(CalcCellFissionSourceNormTest, bart::testing::AllDimensions);

TYPED_TEST(CalcCellFissionSourceNormTest, Constructor) {
  static constexpr int dim = this->dim;

  EXPECT_CALL(*this->finite_element_ptr_, n_cell_quad_pts())
      .WillOnce(DoDefault());

  calculator::cell::FissionSourceNorm<dim> test_calculator(
      this->finite_element_ptr_,
      this->cross_sections_ptr_);

  EXPECT_EQ(this->finite_element_ptr_.use_count(), 2);
  EXPECT_EQ(this->cross_sections_ptr_.use_count(), 2);
}

TYPED_TEST(CalcCellFissionSourceNormTest, GetCellNormNonFissile) {
  static constexpr int dim = this->dim;
  auto& finite_element_mock = *(this->finite_element_ptr_);

  for (auto cell : this->cells_) {
    cell->set_material_id(0);
    EXPECT_CALL(finite_element_mock, SetCell(cell));

    calculator::cell::FissionSourceNorm<dim> test_calculator(
        this->finite_element_ptr_,
        this->cross_sections_ptr_);

    double fission_norm = test_calculator.GetCellNorm(cell);
    EXPECT_EQ(fission_norm, 0.0);
    break;
  }
}



} // namespace