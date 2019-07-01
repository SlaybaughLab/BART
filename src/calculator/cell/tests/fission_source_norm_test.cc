#include "calculator/cell/fission_source_norm.h"

#include <memory>

#include "domain/tests/finite_element_mock.h"
#include "test_helpers/gmock_wrapper.h"


namespace  {

using namespace bart;

template <typename DimensionWrapper>
class CalcCellFissionSourceNormTest :public ::testing::Test {
 protected:
  static constexpr int dim = DimensionWrapper::value;
  using CellPtr = typename domain::FiniteElementI<dim>::CellPtr;
  using FiniteElementType = typename domain::FiniteElementMock<dim>;
  using FissionSourceNormType = typename calculator::cell::FissionSourceNorm<dim>;


  CellPtr cell_ptr_;

  // Supporting mocks
  std::shared_ptr<FiniteElementType> finite_element_ptr_;

  void SetUp() override;
};

template <typename DimensionWrapper>
void CalcCellFissionSourceNormTest<DimensionWrapper>::SetUp() {
  finite_element_ptr_ = std::make_shared<FiniteElementType>();
}

TYPED_TEST_CASE(CalcCellFissionSourceNormTest, bart::testing::AllDimensions);

TYPED_TEST(CalcCellFissionSourceNormTest, Constructor) {
  static constexpr int dim = this->dim;

  calculator::cell::FissionSourceNorm<dim> test_calculator(this->finite_element_ptr_);

  EXPECT_EQ(this->finite_element_ptr_.use_count(), 2);
}

} // namespace