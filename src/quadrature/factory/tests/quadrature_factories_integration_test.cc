#include "quadrature/factory/quadrature_factories.h"

#include "quadrature/ordinate.h"

#include "test_helpers/gmock_wrapper.h"


namespace  {

using namespace bart;

template <typename DimensionWrapper>
class QuadratureFactoriesIntegrationTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
};

TYPED_TEST_CASE(QuadratureFactoriesIntegrationTest, bart::testing::AllDimensions);

/* Call to MakeOrdinate specifying default implementation of OrdinateI should
 * return the correct type.
 */
TYPED_TEST(QuadratureFactoriesIntegrationTest, MakeOrdinate) {
  constexpr int dim = this->dim;
  auto ordinate_ptr = quadrature::factory::MakeOrdinatePtr<dim>(
      quadrature::OrdinateType::kDefault);

  ASSERT_NE(nullptr, ordinate_ptr);

  using ExpectedType = quadrature::Ordinate<dim>;
  EXPECT_NE(nullptr, dynamic_cast<ExpectedType*>(ordinate_ptr.get()));
}


} // namespace