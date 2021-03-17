#include "quadrature/ordinate.hpp"

#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

/* Tests for verifying the operation of the default implementation of the
 * quadrature::OrdianteI class. No mocks are required and tests are performed
 * in all three dimensions. */
template <typename DimensionWrapper>
class QuadratureOrdinateTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
};

TYPED_TEST_CASE(QuadratureOrdinateTest, bart::testing::AllDimensions);

// Constructor should correctly set the cartesian position if provided.
TYPED_TEST(QuadratureOrdinateTest, Construction) {
  constexpr int dim = this->dim;

  std::array<double, dim> position;
  auto random_position = test_helpers::RandomVector(dim, -10, 10);
  for (int i = 0; i < dim; ++i)
    position.at(i) = random_position.at(i);

  quadrature::Ordinate<dim> test_ordinate{quadrature::CartesianPosition<dim>(position)};

  EXPECT_EQ(position, test_ordinate.cartesian_position());

  auto tensor = test_ordinate.cartesian_position_tensor();

  for (int i = 0; i < dim; ++i)
    EXPECT_EQ(position.at(i), tensor[i]);
}
// Default constructor should initialize position as the origin.
TYPED_TEST(QuadratureOrdinateTest, DefaultConstruction) {
  constexpr int dim = this->dim;

  std::array<double, dim> position;
  position.fill(0);

  quadrature::Ordinate<dim> test_ordinate;

  EXPECT_EQ(position, test_ordinate.cartesian_position());
}
// SetCartesianPosition should update stored cartesian position
TYPED_TEST(QuadratureOrdinateTest, SetCartesianPosition) {
  constexpr int dim = this->dim;

  std::array<double, dim> position, new_position;
  auto random_position = test_helpers::RandomVector(dim, -10, 10);
  auto second_position = test_helpers::RandomVector(dim, -10, 10);

  for (int i = 0; i < dim; ++i) {
    position.at(i) = random_position.at(i);
    new_position.at(i) = second_position.at(i);
  }

  quadrature::Ordinate<dim> test_ordinate{quadrature::CartesianPosition<dim>(position)};
  EXPECT_EQ(position, test_ordinate.cartesian_position());

  test_ordinate.set_cartesian_position(quadrature::CartesianPosition<dim>(new_position));
  EXPECT_EQ(new_position, test_ordinate.cartesian_position());
}
// Equality operator should return if cartesian positions are equal
TYPED_TEST(QuadratureOrdinateTest, OperatorEquality) {
  constexpr int dim = this->dim;

  std::array<double, dim> position, negative_position, second_position;
  auto random_position = test_helpers::RandomVector(dim, -10, 10);
  auto random_second_position = test_helpers::RandomVector(dim, -10, 10);
  for (int i = 0; i < dim; ++i) {
    position.at(i) = random_position.at(i);
    negative_position.at(i) = -random_position.at(i);
    second_position.at(i) = random_second_position.at(i);
  }

  quadrature::Ordinate<dim> test_ordinate{quadrature::CartesianPosition<dim>(position)};
  quadrature::Ordinate<dim> test_ordinate_negative{quadrature::CartesianPosition<dim>(negative_position)};
  quadrature::Ordinate<dim> test_ordinate_second{quadrature::CartesianPosition<dim>(second_position)};

  std::array<quadrature::Ordinate<dim>, 3> ordinates{test_ordinate, test_ordinate_negative, test_ordinate_second};

  for (unsigned int i = 0; i < ordinates.size(); ++i) {
    for (unsigned int j = 0; j < ordinates.size(); ++j) {
      EXPECT_EQ(ordinates.at(i) == ordinates.at(j), i == j);
      EXPECT_EQ(ordinates.at(i) == ordinates.at(j).cartesian_position(), i == j);
      EXPECT_EQ(ordinates.at(i) != ordinates.at(j), i != j);
      EXPECT_EQ(ordinates.at(i) != ordinates.at(j).cartesian_position(), i != j);
    }
  }
}

// Reflect free-function should return a reflected ordinate
TYPED_TEST(QuadratureOrdinateTest, Reflect) {
  constexpr int dim = this->dim;

  std::array<double, dim> position, negative_position;
  auto random_position = test_helpers::RandomVector(dim, -10, 10);
  for (int i = 0; i < dim; ++i) {
    position.at(i) = random_position.at(i);
    negative_position.at(i) = -random_position.at(i);
  }

  quadrature::Ordinate<dim> test_ordinate{quadrature::CartesianPosition<dim>(position)};
  quadrature::Ordinate<dim> test_ordinate_negative{quadrature::CartesianPosition<dim>(negative_position)};

  EXPECT_TRUE(test_ordinate == quadrature::Reflect(test_ordinate_negative));
  EXPECT_TRUE(test_ordinate_negative == quadrature::Reflect(test_ordinate));
}


} // namespace