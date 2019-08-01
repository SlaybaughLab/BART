#include "framework/builder/cfem_framework_builder.h"

#include <functional>
#include <numeric>
#include <vector>

#include <deal.II/fe/fe_q.h>

#include "problem/tests/parameters_mock.h"
#include "domain/definition.h"
#include "domain/finite_element_gaussian.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

using ::testing::DoDefault, ::testing::NiceMock, ::testing::Return;

template <typename DimensionWrapper>
class IntegrationTestCFEMFrameworkBuilder : public ::testing::Test {
 protected:
  static constexpr int dim = DimensionWrapper::value;

  using FrameworkBuilder = framework::builder::CFEM_FrameworkBuilder<dim>;
  using ProblemParameters = NiceMock<problem::ParametersMock>;

  FrameworkBuilder test_builder;
  ProblemParameters parameters;

  // Test Parameters
  const int polynomial_degree = 2;
  std::vector<double> spatial_max;
  std::vector<int> n_cells;
  std::array<int, 4> dofs_per_cell_by_dim_{1, 3, 9, 27};

  void SetUp() override;

};

template <typename DimensionWrapper>
void IntegrationTestCFEMFrameworkBuilder<DimensionWrapper>::SetUp() {
  for (int i = 0; i < this->dim; ++i) {
    spatial_max.push_back(i + 10);
    n_cells.push_back(i + 10);
  }

  ON_CALL(parameters, NCells())
      .WillByDefault(Return(n_cells));
  ON_CALL(parameters, SpatialMax())
      .WillByDefault(Return(spatial_max));
  ON_CALL(parameters, FEPolynomialDegree())
      .WillByDefault(Return(polynomial_degree));

}

TYPED_TEST_CASE(IntegrationTestCFEMFrameworkBuilder,
                bart::testing::AllDimensions);

TYPED_TEST(IntegrationTestCFEMFrameworkBuilder, BuildFiniteElementTest) {
  EXPECT_CALL(this->parameters, FEPolynomialDegree())
      .WillOnce(DoDefault());

  auto finite_element_ptr =
      this->test_builder.BuildFiniteElement(&this->parameters);

  // Test that correct type was returned
  auto gaussian_ptr = dynamic_cast<domain::FiniteElementGaussian<this->dim>*>(
      finite_element_ptr.get());
  EXPECT_NE(gaussian_ptr, nullptr);

  // Test that correct values were used in instantiation
  EXPECT_EQ(finite_element_ptr->polynomial_degree(), this->polynomial_degree);
  auto dealii_finite_element_ptr = dynamic_cast<dealii::FE_Q<this->dim>*>(
      finite_element_ptr->finite_element());
  EXPECT_NE(dealii_finite_element_ptr, nullptr);
}

TYPED_TEST(IntegrationTestCFEMFrameworkBuilder, BuildDomainTest) {

  auto finite_element_ptr =
      this->test_builder.BuildFiniteElement(&this->parameters);

  EXPECT_CALL(this->parameters, NCells())
      .WillOnce(DoDefault());
  EXPECT_CALL(this->parameters, SpatialMax())
      .WillOnce(DoDefault());

  auto test_domain_ptr =
      this->test_builder.BuildDomain(
          &this->parameters, finite_element_ptr, "1 1 2 2");

  EXPECT_NE(dynamic_cast<domain::Definition<this->dim>*>(test_domain_ptr.get()),
      nullptr);
}

} // namespace