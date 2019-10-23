#include "quadrature/utility/quadrature_utilities.h"

#include "quadrature/tests/quadrature_generator_mock.h"
#include "quadrature/tests/quadrature_set_mock.h"
#include "quadrature/tests/ordinate_mock.h"

#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class QuadratureUtilityTests : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
};

TYPED_TEST_CASE(QuadratureUtilityTests, bart::testing::AllDimensions);

TYPED_TEST(QuadratureUtilityTests, Dummy) {
  EXPECT_TRUE(true);
}

TYPED_TEST(QuadratureUtilityTests, Reflect) {
  constexpr int dim = this->dim;

  std::array<double, dim> position, negative_position;
  auto random_position = btest::RandomVector(dim, -10, 10);
  for (int i = 0; i < dim; ++i) {
    position.at(i) = random_position.at(i);
    negative_position.at(i) = -random_position.at(i);
  }

  quadrature::OrdinateMock<dim> mock_ordinate;

  EXPECT_CALL(mock_ordinate, cartesian_position())
      .WillOnce(::testing::Return(position));

  EXPECT_EQ(quadrature::utility::Reflect(mock_ordinate), negative_position);

  EXPECT_CALL(mock_ordinate, cartesian_position())
      .WillOnce(::testing::Return(negative_position));

  EXPECT_EQ(quadrature::utility::Reflect(mock_ordinate), position);
}


} // namespace