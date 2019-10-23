#include "quadrature/utility/quadrature_utilities.h"

#include <algorithm>
#include <cmath>

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

TYPED_TEST(QuadratureUtilityTests, ReflectAcrossOrigin) {
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

  EXPECT_EQ(quadrature::utility::ReflectAcrossOrigin(mock_ordinate),
            negative_position);

  EXPECT_CALL(mock_ordinate, cartesian_position())
      .WillOnce(::testing::Return(negative_position));

  EXPECT_EQ(quadrature::utility::ReflectAcrossOrigin(mock_ordinate),
            position);
}

TYPED_TEST(QuadratureUtilityTests, GenerateAllPositiveX) {
  const int dim = this->dim;
  const int n_points = 3;
  const int n_quadrants = std::pow(2, dim)/2;

  // Quadrature set to be distributed in positive x
  std::vector<std::pair<quadrature::CartesianPosition<dim>, quadrature::Weight>>
      quadrature_set;

  for (int i = 0; i < n_points; ++i) {
    auto random_position = btest::RandomVector(dim, 1, 10);
    auto random_weight = btest::RandomDouble(0, 2);
    std::array<double, dim> position;
    for (int j = 0; j < dim; ++j)
      position.at(j) = random_position.at(j);
    quadrature_set.emplace_back(quadrature::CartesianPosition<dim>(position),
                                quadrature::Weight(random_weight));
  }

  auto distributed_set = quadrature::utility::GenerateAllPositiveX<dim>(quadrature_set);

  EXPECT_EQ(distributed_set.size(), n_points * n_quadrants);

  for (const auto& quadrature_pair : quadrature_set) {
    EXPECT_EQ(1, std::count(distributed_set.cbegin(), distributed_set.cend(),
                            quadrature_pair));
    if (dim > 1) {
      auto negative_y_pair = quadrature_pair;
      negative_y_pair.first.get().at(1) *= -1;
      EXPECT_EQ(1, std::count(distributed_set.cbegin(), distributed_set.cend(),
                              negative_y_pair));
    }
    if (dim > 2) {
      auto negative_z_pair = quadrature_pair;

      negative_z_pair.first.get().at(2) *= -1;
      EXPECT_EQ(1, std::count(distributed_set.cbegin(), distributed_set.cend(),
                              negative_z_pair));
      auto negative_zy_pair = negative_z_pair;
      negative_zy_pair.first.get().at(1) *= -1;
      EXPECT_EQ(1, std::count(distributed_set.cbegin(), distributed_set.cend(),
                              negative_zy_pair));

    }
  }
}




} // namespace