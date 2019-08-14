#include "domain/finite_element/finite_element_factory.h"

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>

#include "domain/finite_element/finite_element_gaussian.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class DomainFiniteElementFactoryIntegrationTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
};

TYPED_TEST_CASE(DomainFiniteElementFactoryIntegrationTest,
                bart::testing::AllDimensions);

TYPED_TEST(DomainFiniteElementFactoryIntegrationTest,
    MakeContinuousFiniteElement) {
  constexpr int dim = this->dim;
  const auto discretization_type = problem::DiscretizationType::kContinuousFEM;
  const int polynomial_degree = 2;
  const auto implementation_type =
      domain::finite_element::FiniteElementImpl::kGaussian;

  auto finite_element_ptr =
      domain::finite_element::FiniteElementFactory<dim>::MakeFiniteElement(
          discretization_type, polynomial_degree, implementation_type);

  ASSERT_NE(finite_element_ptr, nullptr);

  using ExpectedType = domain::finite_element::FiniteElementGaussian<dim>;
  ASSERT_NE(nullptr, dynamic_cast<ExpectedType*>(finite_element_ptr.get()));
  EXPECT_EQ(finite_element_ptr->polynomial_degree(), polynomial_degree);

  using ExpectedDealiiType = dealii::FE_Q<dim>;

  EXPECT_NE(nullptr,
            dynamic_cast<ExpectedDealiiType*>(
                finite_element_ptr->finite_element()));
}

TYPED_TEST(DomainFiniteElementFactoryIntegrationTest,
           MakeDiscontinuousFiniteElement) {
  constexpr int dim = this->dim;
  const auto discretization_type = problem::DiscretizationType::kDiscontinuousFEM;
  const int polynomial_degree = 2;
  const auto implementation_type =
      domain::finite_element::FiniteElementImpl::kGaussian;

  auto finite_element_ptr =
      domain::finite_element::FiniteElementFactory<dim>::MakeFiniteElement(
          discretization_type, polynomial_degree, implementation_type);

  ASSERT_NE(finite_element_ptr, nullptr);

  using ExpectedType = domain::finite_element::FiniteElementGaussian<dim>;
  ASSERT_NE(nullptr, dynamic_cast<ExpectedType*>(finite_element_ptr.get()));
  EXPECT_EQ(finite_element_ptr->polynomial_degree(), polynomial_degree);

  using ExpectedDealiiType = dealii::FE_DGQ<dim>;

  EXPECT_NE(nullptr,
            dynamic_cast<ExpectedDealiiType*>(
                finite_element_ptr->finite_element()));
}

} // namespace