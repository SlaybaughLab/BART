#include "calculator/cell/fission_source_norm.h"

#include <memory>

#include "data/cross_sections.h"
#include "domain/tests/finite_element_mock.h"
#include "material/tests/mock_material.h"
#include "test_helpers/gmock_wrapper.h"


namespace  {

using namespace bart;

using ::testing::NiceMock;

template <typename DimensionWrapper>
class CalcCellFissionSourceNormTest :public ::testing::Test {
 protected:
  static constexpr int dim = DimensionWrapper::value;
  using CellPtr = typename domain::FiniteElementI<dim>::CellPtr;
  using FiniteElementType = typename domain::FiniteElementMock<dim>;
  using FissionSourceNormType = typename calculator::cell::FissionSourceNorm<dim>;

  CellPtr cell_ptr_;

  // Supporting objects and mocks
  std::shared_ptr<FiniteElementType> finite_element_ptr_;
  std::shared_ptr<data::CrossSections> cross_sections_ptr_;


  NiceMock<btest::MockMaterial> mock_material_;

  void SetUp() override;
};

template <typename DimensionWrapper>
void CalcCellFissionSourceNormTest<DimensionWrapper>::SetUp() {
  finite_element_ptr_ = std::make_shared<FiniteElementType>();
  cross_sections_ptr_ = std::make_shared<data::CrossSections>(mock_material_);
}

TYPED_TEST_CASE(CalcCellFissionSourceNormTest, bart::testing::AllDimensions);

TYPED_TEST(CalcCellFissionSourceNormTest, Constructor) {
  static constexpr int dim = this->dim;

  calculator::cell::FissionSourceNorm<dim> test_calculator(
      this->finite_element_ptr_,
      this->cross_sections_ptr_);

  EXPECT_EQ(this->finite_element_ptr_.use_count(), 2);
  EXPECT_EQ(this->cross_sections_ptr_.use_count(), 2);
}

} // namespace