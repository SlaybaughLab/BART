#include "calculator/residual/tests/cell_isotropic_residual_mock.hpp"
#include "calculator/residual/domain_isotropic_residual.hpp"
#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_assertions.hpp"

namespace  {

using namespace bart;

using ::testing::DoDefault, ::testing::NiceMock, ::testing::Return, ::testing::_;

template <typename DimensionWrapper>
class CalculatorResidualDomainIsotropicResidualTest : public ::testing::Test,
                                                      public bart::testing::DealiiTestDomain<DimensionWrapper::value> {
 public:
  static constexpr int dim{ DimensionWrapper::value };
};

TYPED_TEST_SUITE(CalculatorResidualDomainIsotropicResidualTest, bart::testing::AllDimensions);

TYPED_TEST(CalculatorResidualDomainIsotropicResidualTest, Dummy) {
  EXPECT_TRUE(false);
}

} // namespace
