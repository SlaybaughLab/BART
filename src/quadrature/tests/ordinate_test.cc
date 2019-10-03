#include "quadrature/ordinate.h"

#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class QuadratureOrdinateTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
};


TYPED_TEST_CASE(QuadratureOrdinateTest, bart::testing::AllDimensions);

TYPED_TEST(QuadratureOrdinateTest, Construction) {
  constexpr int dim = this->dim;

  std::array<double, dim> position;
  auto random_position = btest::RandomVector(dim, -10, 10);
  for (int i = 0; i < dim; ++i)
    position.at(i) = random_position.at(i);

  quadrature::Ordinate<dim> test_ordinate{quadrature::CartesianPosition<dim>(position)};

  EXPECT_EQ(position, test_ordinate.cartesian_position());
}


} // namespace