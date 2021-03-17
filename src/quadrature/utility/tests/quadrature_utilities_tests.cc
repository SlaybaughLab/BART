#include "quadrature/utility/quadrature_utilities.h"

#include <algorithm>
#include <cmath>

#include "quadrature/tests/quadrature_point_mock.hpp"
#include "quadrature/tests/quadrature_generator_mock.h"
#include "quadrature/tests/quadrature_set_mock.h"
#include "quadrature/tests/ordinate_mock.hpp"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

using namespace bart;

/* Tests for various quadrature utilities. Will be run in all three dimensions.
 */
template <typename DimensionWrapper>
class QuadratureUtilityTests : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
};

TYPED_TEST_CASE(QuadratureUtilityTests, bart::testing::AllDimensions);

// ReflectAcrossOrigin should return the correct reflected value
TYPED_TEST(QuadratureUtilityTests, ReflectAcrossOrigin) {
  constexpr int dim = this->dim;
  std::array<double, dim> position, negative_position;
  auto random_position = test_helpers::RandomVector(dim, -10, 10);

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

// GenerateAllPositiveX should generate all correct points
TYPED_TEST(QuadratureUtilityTests, GenerateAllPositiveX) {
  const int dim = this->dim;
  const int n_points = 3;
  const int n_quadrants = std::pow(2, dim)/2;

  std::vector<std::pair<quadrature::CartesianPosition<dim>, quadrature::Weight>>
      quadrature_set;

  // Fill the quadrature set with random positions
  for (int i = 0; i < n_points; ++i) {
    auto random_position = test_helpers::RandomVector(dim, 1, 10);
    auto random_weight = test_helpers::RandomDouble(0, 2);
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

/* GenerateAllPositiveXScalar shouldn't have any effect on quadrature points
 * at the origin (the only points without a reflection) */
TYPED_TEST(QuadratureUtilityTests, GenerateAllPositiveXScalar) {
  const int dim = this->dim;

  std::vector<std::pair<quadrature::CartesianPosition<dim>, quadrature::Weight>>
      quadrature_set;
  std::array<double, dim> zero_position;
  zero_position.fill(0);
  quadrature_set.emplace_back(quadrature::CartesianPosition<dim>(zero_position),
                              quadrature::Weight(1));

  auto distributed_set =
      quadrature::utility::GenerateAllPositiveX<dim>(quadrature_set);

  EXPECT_EQ(distributed_set.size(), 1);
  EXPECT_EQ(quadrature_set.at(0), distributed_set.at(0));
}

// QuadraturePointCompare should provide a correct less-than operation.
TYPED_TEST(QuadratureUtilityTests, QuadraturePointCompare) {
  const int dim = this->dim;

  std::array<double, dim> position_1, position_2;
  auto random_position_1 = test_helpers::RandomVector(dim, -100, 100);
  auto random_position_2 = test_helpers::RandomVector(dim, -100, 100);

  for (int i = 0; i < dim; ++i) {
    position_1.at(i) = random_position_1.at(i);
    position_2.at(i) = random_position_2.at(i);
  }

  bool expected_result = (position_1 < position_2);
  auto mock_point_1 = std::make_shared<quadrature::QuadraturePointMock<dim>>();
  auto mock_point_2 = std::make_shared<quadrature::QuadraturePointMock<dim>>();

  EXPECT_CALL(*mock_point_1, cartesian_position())
      .WillOnce(::testing::Return(position_1));
  EXPECT_CALL(*mock_point_2, cartesian_position())
      .WillOnce(::testing::Return(position_2));

  quadrature::utility::quadrature_point_compare<dim> compare_struct;
  auto result = compare_struct(mock_point_1, mock_point_2);
  EXPECT_EQ(result, expected_result);

  // Check const operations
  const auto& const_compare_struct = compare_struct;
  EXPECT_CALL(*mock_point_1, cartesian_position())
      .WillOnce(::testing::Return(position_1));
  EXPECT_CALL(*mock_point_2, cartesian_position())
      .WillOnce(::testing::Return(position_2));

  auto const_result = const_compare_struct(mock_point_1, mock_point_2);
  EXPECT_EQ(const_result, expected_result);
}

} // namespace