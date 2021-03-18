#include "quadrature/quadrature_set.hpp"

#include <memory>

#include "quadrature/tests/quadrature_point_mock.hpp"
#include "quadrature/quadrature_types.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

using ::testing::NiceMock;
using ::testing::Return;
/* Tests operation of the default implementation of the quadrature set. This is
 * accomplished using mock quadrature points. Test is performed in all three
 * dimensions.
 *
 * SetUp: Creates mock quadrature points and adds two of them to the quadrature
 * set.
 *
 */
template <typename DimensionWrapper>
class QuadratureSetTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using QuadraturePointMock = NiceMock<quadrature::QuadraturePointMock<dim>>;
  // Object to test
  quadrature::QuadratureSet<dim> test_set_;
  // Pointers to hold mock quadrature points
  std::shared_ptr<QuadraturePointMock> quadrature_point_{ std::make_shared<QuadraturePointMock>()};
  std::shared_ptr<QuadraturePointMock> second_quadrature_point_{ std::make_shared<QuadraturePointMock>()};
  std::shared_ptr<QuadraturePointMock> third_quadrature_point_{ std::make_shared<QuadraturePointMock>()};

  /* Sets up mock quadrature points to return correct positions and add
   * the first and second mock point to the set.
   */
  auto SetUp() -> void override;
};

template <typename DimensionWrapper>
auto QuadratureSetTest<DimensionWrapper>::SetUp() -> void {
  std::array<double, dim> position_1, position_2, position_3;
  position_1.fill(1); position_2.fill(2); position_3.fill(3);

  ON_CALL(*quadrature_point_, cartesian_position()).WillByDefault(::testing::Return(position_1));
  ON_CALL(*second_quadrature_point_, cartesian_position()).WillByDefault(::testing::Return(position_2));
  ON_CALL(*third_quadrature_point_, cartesian_position()).WillByDefault(::testing::Return(position_3));

  test_set_.AddPoint(quadrature_point_);
  test_set_.AddPoint(second_quadrature_point_);
}

TYPED_TEST_CASE(QuadratureSetTest, bart::testing::AllDimensions);

// Default constructor should not throw.
TYPED_TEST(QuadratureSetTest, Constructor) {
  EXPECT_NO_THROW(quadrature::QuadratureSet<this->dim> test_set);
}

// Add point should properly insert and increase the size of the quadrature set.
TYPED_TEST(QuadratureSetTest, AddPoint) {
  constexpr int dim = this->dim;
  quadrature::QuadratureSet<dim> test_set;
  EXPECT_EQ(test_set.size(), 0);
  // Insert one point
  EXPECT_TRUE(test_set.AddPoint(this->quadrature_point_));
  EXPECT_EQ(test_set.size(), 1);
  // Insert second point
  EXPECT_TRUE(test_set.AddPoint(this->second_quadrature_point_));
  EXPECT_EQ(test_set.size(), 2);
  // Attempt insertion of first point again
  EXPECT_FALSE(test_set.AddPoint(this->quadrature_point_));
  EXPECT_EQ(test_set.size(), 2);
}

// Add point should have the appropriate entry in the quadrature_indices_
TYPED_TEST(QuadratureSetTest, AddPointIndices) {
  constexpr int dim = this->dim;
  quadrature::QuadratureSet<dim> test_set;
  std::set<int> expected_indices = {};
  EXPECT_EQ(test_set.quadrature_point_indices().size(), 0);
  // Insert one point
  test_set.AddPoint(this->quadrature_point_);
  expected_indices.insert(0);
  EXPECT_THAT(test_set.quadrature_point_indices(), ::testing::ContainerEq(expected_indices));
  // Insert second point
  test_set.AddPoint(this->second_quadrature_point_);
  expected_indices.insert(1);
  EXPECT_THAT(test_set.quadrature_point_indices(), ::testing::ContainerEq(expected_indices));
  // Attempt insertion of first point again
  test_set.AddPoint(this->quadrature_point_);
  EXPECT_THAT(test_set.quadrature_point_indices(), ::testing::ContainerEq(expected_indices));
}

TYPED_TEST(QuadratureSetTest, GetBoundaryReflectionBadReflection) {
  using Boundary = problem::Boundary;
  EXPECT_ANY_THROW({
    [[maybe_unused]] auto boundary = this->test_set_.GetBoundaryReflection(this->quadrature_point_, Boundary::kXMin);
  });
}

TYPED_TEST(QuadratureSetTest, GetBoundaryReflection) {
  constexpr int dim = this->dim;
  using MockQuadraturePointType =  NiceMock<quadrature::QuadraturePointMock<dim>>;
  using Boundary = problem::Boundary;
  // Create reflected quadrature points and add them
  auto mock_x_reflection = std::make_shared<MockQuadraturePointType>();

  std::array<double, dim> x_reflected_position;
  x_reflected_position.fill(1);
  x_reflected_position.at(0) = -1;
  ON_CALL(*mock_x_reflection, cartesian_position()).WillByDefault(Return(x_reflected_position));

  this->test_set_.AddPoint(mock_x_reflection);
  EXPECT_EQ(this->test_set_.GetBoundaryReflection(this->quadrature_point_, Boundary::kXMax), mock_x_reflection);
  EXPECT_EQ(this->test_set_.GetBoundaryReflection(mock_x_reflection, Boundary::kXMin), this->quadrature_point_);

  if (dim > 1) {
    std::array<double, dim> y_reflected_position;
    y_reflected_position.fill(1);
    y_reflected_position.at(1) = -1;
    auto mock_y_reflection = std::make_shared<MockQuadraturePointType>();
    ON_CALL(*mock_y_reflection, cartesian_position()).WillByDefault(Return(y_reflected_position));
    this->test_set_.AddPoint(mock_y_reflection);
    EXPECT_EQ(this->test_set_.GetBoundaryReflection(this->quadrature_point_, Boundary::kYMax), mock_y_reflection);
    EXPECT_EQ(this->test_set_.GetBoundaryReflection(mock_y_reflection, Boundary::kYMin), this->quadrature_point_);
  }
  if (dim > 2) {
    std::array<double, dim> z_reflected_position;
    z_reflected_position.fill(1);
    z_reflected_position.at(2) = -1;
    auto mock_z_reflection = std::make_shared<MockQuadraturePointType>();
    ON_CALL(*mock_z_reflection, cartesian_position()).WillByDefault(Return(z_reflected_position));
    this->test_set_.AddPoint(mock_z_reflection);
    EXPECT_EQ(this->test_set_.GetBoundaryReflection(this->quadrature_point_, Boundary::kZMax), mock_z_reflection);
    EXPECT_EQ(this->test_set_.GetBoundaryReflection(mock_z_reflection, Boundary::kZMin), this->quadrature_point_);
  }
}

// Getters for quadrature point and index should retrieve the correct values
TYPED_TEST(QuadratureSetTest, GetQuadraturePointAndIndex) {
  constexpr int dim = this->dim;
  quadrature::QuadratureSet<dim> test_set;
  test_set.AddPoint(this->quadrature_point_);
  // Add first point
  EXPECT_EQ(this->quadrature_point_, test_set.GetQuadraturePoint(quadrature::QuadraturePointIndex(0)));
  EXPECT_EQ(test_set.GetQuadraturePointIndex(this->quadrature_point_), 0);
  // Add second point
  test_set.AddPoint(this->second_quadrature_point_);
  EXPECT_EQ(this->second_quadrature_point_, test_set.GetQuadraturePoint(quadrature::QuadraturePointIndex(1)));
  EXPECT_EQ(test_set.GetQuadraturePointIndex(this->second_quadrature_point_), 1);
}

// Trying to add a null quadrature point ptr should throw
TYPED_TEST(QuadratureSetTest, AddPointNullPtr) {
  quadrature::QuadratureSet<this->dim> test_set;
  EXPECT_ANY_THROW(test_set.AddPoint(nullptr));
}

/* Points in set without reflection will return nullptr when the reflection is
 * requested and index will return an empty optional.
 */
TYPED_TEST(QuadratureSetTest, DefaultGetReflection) {
  EXPECT_EQ(nullptr, this->test_set_.GetReflection(this->quadrature_point_));
  EXPECT_FALSE(this->test_set_.GetReflectionIndex(this->quadrature_point_).has_value());
  EXPECT_EQ(nullptr, this->test_set_.GetReflection(this->second_quadrature_point_));
  EXPECT_FALSE(this->test_set_.GetReflectionIndex(this->second_quadrature_point_).has_value());
  EXPECT_ANY_THROW([[maybe_unused]] auto reflection = this->test_set_.GetReflection(this->third_quadrature_point_));
  EXPECT_ANY_THROW(this->test_set_.GetReflectionIndex(this->third_quadrature_point_).has_value());
}

// SetReflection should set the reflection properly.
TYPED_TEST(QuadratureSetTest, SetReflection) {
  this->test_set_.SetReflection(this->quadrature_point_, this->second_quadrature_point_);
  EXPECT_EQ(this->test_set_.GetReflection(this->quadrature_point_), this->second_quadrature_point_);
  EXPECT_EQ(this->test_set_.GetReflection(this->second_quadrature_point_), this->quadrature_point_);

  EXPECT_EQ(this->test_set_.GetReflectionIndex(this->quadrature_point_).value_or(-1),
            this->test_set_.GetQuadraturePointIndex(this->second_quadrature_point_));
  EXPECT_EQ(this->test_set_.GetReflectionIndex(this->second_quadrature_point_).value_or(-1),
            this->test_set_.GetQuadraturePointIndex(this->quadrature_point_));

  // Set again and reverse, verify nothing has changed, returns false because no
  // insertion.
  this->test_set_.SetReflection(this->quadrature_point_, this->second_quadrature_point_);
  this->test_set_.SetReflection(this->second_quadrature_point_, this->quadrature_point_);
  EXPECT_EQ(this->test_set_.GetReflection(this->quadrature_point_), this->second_quadrature_point_);
  EXPECT_EQ(this->test_set_.GetReflection(this->second_quadrature_point_), this->quadrature_point_);
  EXPECT_EQ(this->test_set_.GetReflectionIndex(this->quadrature_point_).value_or(-1),
            this->test_set_.GetQuadraturePointIndex(this->second_quadrature_point_));
  EXPECT_EQ(this->test_set_.GetReflectionIndex(this->second_quadrature_point_).value_or(-1),
            this->test_set_.GetQuadraturePointIndex(this->quadrature_point_));
}

// SetReflection should throw if a point is being set as its own reflection.
TYPED_TEST(QuadratureSetTest, SamePointAsReflection) {
  EXPECT_ANY_THROW({
    this->test_set_.SetReflection(this->quadrature_point_, this->quadrature_point_);
  });
}

// SetReflection cannot be called
TYPED_TEST(QuadratureSetTest, ReflectionThirdPoint) {
  EXPECT_ANY_THROW({ this->test_set_.SetReflection(this->third_quadrature_point_, this->quadrature_point_); });
  EXPECT_ANY_THROW({ this->test_set_.SetReflection(this->quadrature_point_, this->third_quadrature_point_); });
  EXPECT_ANY_THROW({ this->test_set_.SetReflection(this->third_quadrature_point_, this->third_quadrature_point_); });
}

/* If you replace a reflection of a point, the original reflection should no
 * longer consider the first point to be its reflection.
 */
TYPED_TEST(QuadratureSetTest, SetNewReflectionFirstPoint) {
  this->test_set_.AddPoint(this->third_quadrature_point_);
  this->test_set_.SetReflection(this->quadrature_point_, this->second_quadrature_point_);
  this->test_set_.SetReflection(this->quadrature_point_, this->third_quadrature_point_);
  // Second point now has no reflection
  EXPECT_EQ(nullptr, this->test_set_.GetReflection(this->second_quadrature_point_));
  EXPECT_FALSE(this->test_set_.GetReflectionIndex(this->second_quadrature_point_).has_value());
  // Third and first now reflect each other
  EXPECT_EQ(this->test_set_.GetReflection(this->quadrature_point_), this->third_quadrature_point_);
  EXPECT_EQ(this->test_set_.GetReflection(this->third_quadrature_point_), this->quadrature_point_);
  EXPECT_EQ(this->test_set_.GetReflectionIndex(this->quadrature_point_).value_or(-1),
            this->test_set_.GetQuadraturePointIndex(this->third_quadrature_point_));
  EXPECT_EQ(this->test_set_.GetReflectionIndex(this->third_quadrature_point_).value_or(-1),
            this->test_set_.GetQuadraturePointIndex(this->quadrature_point_));
}

/* Same as previous test but with second point */
TYPED_TEST(QuadratureSetTest, SetNewReflectionSecondPoint) {
  this->test_set_.AddPoint(this->third_quadrature_point_);

  this->test_set_.SetReflection(this->quadrature_point_, this->second_quadrature_point_);
  this->test_set_.SetReflection(this->third_quadrature_point_, this->quadrature_point_);
  // Second point now has no reflection
  EXPECT_EQ(nullptr, this->test_set_.GetReflection(this->second_quadrature_point_));
  EXPECT_FALSE(this->test_set_.GetReflectionIndex(this->second_quadrature_point_).has_value());
  // Third and first now reflect each other
  EXPECT_EQ(this->test_set_.GetReflection(this->quadrature_point_), this->third_quadrature_point_);
  EXPECT_EQ(this->test_set_.GetReflection(this->third_quadrature_point_), this->quadrature_point_);
  EXPECT_EQ(this->test_set_.GetReflectionIndex(this->quadrature_point_).value_or(-1),
            this->test_set_.GetQuadraturePointIndex(this->third_quadrature_point_));
  EXPECT_EQ(this->test_set_.GetReflectionIndex(this->third_quadrature_point_).value_or(-1),
            this->test_set_.GetQuadraturePointIndex(this->quadrature_point_));
}

// Iterator should properly return an interator to the set
TYPED_TEST(QuadratureSetTest, Iterator) {
  ASSERT_EQ(this->test_set_.size(), 2);
  for (auto& point_ptr : this->test_set_) {
    if (point_ptr != this->quadrature_point_){
      EXPECT_EQ(point_ptr, this->second_quadrature_point_);
    } else {
      EXPECT_EQ(point_ptr, this->quadrature_point_);
    }
  }
}

// ConstIterator should properly return an interator to the set
TYPED_TEST(QuadratureSetTest, ConstIterator) {
  ASSERT_EQ(this->test_set_.size(), 2);
  for (auto it = this->test_set_.cbegin(); it != this->test_set_.cend(); ++it) {
    auto point_ptr = *it;
    if (point_ptr != this->quadrature_point_) {
      EXPECT_EQ(point_ptr, this->second_quadrature_point_);
    } else {
      EXPECT_EQ(point_ptr, this->quadrature_point_);
    }
  }
}

} // namespace