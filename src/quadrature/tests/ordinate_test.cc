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

  auto tensor = test_ordinate.cartesian_position_tensor();

  for (int i = 0; i < dim; ++i)
    EXPECT_EQ(position.at(i), tensor[i]);
}

TYPED_TEST(QuadratureOrdinateTest, OperatorEqual) {
  constexpr int dim = this->dim;

  std::array<double, dim> position, negative_position, second_position;
  auto random_position = btest::RandomVector(dim, -10, 10);
  auto random_second_position = btest::RandomVector(dim, -10, 10);
  for (int i = 0; i < dim; ++i) {
    position.at(i) = random_position.at(i);
    negative_position.at(i) = -random_position.at(i);
    second_position.at(i) = random_second_position.at(i);
  }

  quadrature::Ordinate<dim> test_ordinate{
    quadrature::CartesianPosition<dim>(position)};
  quadrature::Ordinate<dim> test_ordinate_negative{
    quadrature::CartesianPosition<dim>(negative_position)};
  quadrature::Ordinate<dim> test_ordinate_second{
      quadrature::CartesianPosition<dim>(second_position)};

  std::array<quadrature::Ordinate<dim>, 3> ordinates{test_ordinate,
                                                     test_ordinate_negative,
                                                     test_ordinate_second};

  for (unsigned int i = 0; i < ordinates.size(); ++i) {
    for (unsigned int j = 0; j < ordinates.size(); ++j) {
      if (i == j) {
        EXPECT_TRUE(ordinates.at(i) == ordinates.at(j))
                  << "i,j = " << i << ", " << j;
        EXPECT_FALSE(ordinates.at(i) != ordinates.at(j))
                  << "i,j = " << i << ", " << j;
      } else {
        EXPECT_TRUE(ordinates.at(i) != ordinates.at(j))
                  << "i,j = " << i << ", " << j;
        EXPECT_FALSE(ordinates.at(i) == ordinates.at(j))
                  << "i,j = " << i << ", " << j;
      }
    }
  }
}


} // namespace