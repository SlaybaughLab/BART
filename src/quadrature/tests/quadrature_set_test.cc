#include "quadrature/quadrature_set.h"

#include <memory>

#include "quadrature/tests/quadrature_point_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

using ::testing::NiceMock;

template <typename DimensionWrapper>
class QuadratureSetTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;

  quadrature::QuadratureSet<dim> test_set_;
  std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point_,
      second_quadrature_point_, third_quadrature_point_;

  void SetUp() override {

    std::array<double, dim> position_1, position_2, position_3;
    position_1.fill(1); position_2.fill(2); position_3.fill(3);

    auto mock_quadrature_point_ =
        std::make_shared<NiceMock<quadrature::QuadraturePointMock<dim>>>();
    auto mock_second_quadrature_point_ =
        std::make_shared<NiceMock<quadrature::QuadraturePointMock<dim>>>();
    auto mock_third_quadrature_point_ =
        std::make_shared<NiceMock<quadrature::QuadraturePointMock<dim>>>();

    ON_CALL(*mock_quadrature_point_, cartesian_position()).WillByDefault(::testing::Return(position_1));
    ON_CALL(*mock_second_quadrature_point_, cartesian_position()).WillByDefault(::testing::Return(position_2));
    ON_CALL(*mock_third_quadrature_point_, cartesian_position()).WillByDefault(::testing::Return(position_3));

    quadrature_point_ = mock_quadrature_point_;
    second_quadrature_point_ = mock_second_quadrature_point_;
    third_quadrature_point_ = mock_third_quadrature_point_;

    test_set_.AddPoint(quadrature_point_);
    test_set_.AddPoint(second_quadrature_point_);
  }
};

TYPED_TEST_CASE(QuadratureSetTest, bart::testing::AllDimensions);

TYPED_TEST(QuadratureSetTest, Constructor) {
  constexpr int dim = this->dim;
  EXPECT_NO_THROW(quadrature::QuadratureSet<dim> test_set);
}

TYPED_TEST(QuadratureSetTest, AddPoint) {
  constexpr int dim = this->dim;
  quadrature::QuadratureSet<dim> test_set;

  EXPECT_EQ(test_set.size(), 0);

  EXPECT_TRUE(test_set.AddPoint(this->quadrature_point_));
  EXPECT_EQ(test_set.size(), 1);

  EXPECT_TRUE(test_set.AddPoint(this->second_quadrature_point_));
  EXPECT_EQ(test_set.size(), 2);

  EXPECT_FALSE(test_set.AddPoint(this->quadrature_point_));
  EXPECT_EQ(test_set.size(), 2);
}

TYPED_TEST(QuadratureSetTest, AddPointNullPtr) {
  constexpr int dim = this->dim;
  quadrature::QuadratureSet<dim> test_set;

  EXPECT_ANY_THROW(test_set.AddPoint(nullptr));
}

TYPED_TEST(QuadratureSetTest, DefaultGetReflection) {
  // Points in set without reflection will return nullptr
  EXPECT_EQ(nullptr, this->test_set_.GetReflection(this->quadrature_point_));
  EXPECT_EQ(nullptr, this->test_set_.GetReflection(this->second_quadrature_point_));
  // Points not in set will throw error
  EXPECT_ANY_THROW(this->test_set_.GetReflection(this->third_quadrature_point_));
}

TYPED_TEST(QuadratureSetTest, SetReflection) {
  // Set and verify
  this->test_set_.SetReflection(this->quadrature_point_,
                                this->second_quadrature_point_);
  EXPECT_EQ(this->test_set_.GetReflection(this->quadrature_point_),
            this->second_quadrature_point_);
  EXPECT_EQ(this->test_set_.GetReflection(this->second_quadrature_point_),
            this->quadrature_point_);

  // Set again and reverse, verify nothing has changed, returns false because no
  // insertion
  this->test_set_.SetReflection(this->quadrature_point_,
                                this->second_quadrature_point_);
  this->test_set_.SetReflection(this->second_quadrature_point_,
                                this->quadrature_point_);
  EXPECT_EQ(this->test_set_.GetReflection(this->quadrature_point_),
            this->second_quadrature_point_);
  EXPECT_EQ(this->test_set_.GetReflection(this->second_quadrature_point_),
            this->quadrature_point_);
}

TYPED_TEST(QuadratureSetTest, SamePointAsReflection) {
  // Same point cannot be its own reflection
  EXPECT_ANY_THROW({
    this->test_set_.SetReflection(this->quadrature_point_,
                                  this->quadrature_point_);
  });
}

TYPED_TEST(QuadratureSetTest, ReflectionThirdPoint) {
  // Point not in set can't be set as a reflection
  EXPECT_ANY_THROW({
    this->test_set_.SetReflection(this->third_quadrature_point_,
                                  this->quadrature_point_);
                   });
  EXPECT_ANY_THROW({
    this->test_set_.SetReflection(this->quadrature_point_,
                                  this->third_quadrature_point_);
                   });
  EXPECT_ANY_THROW({
    this->test_set_.SetReflection(this->third_quadrature_point_,
                                  this->third_quadrature_point_);
                   });
}

TYPED_TEST(QuadratureSetTest, SetNewReflectionFirstPoint) {
  this->test_set_.AddPoint(this->third_quadrature_point_);

  this->test_set_.SetReflection(this->quadrature_point_,
                                this->second_quadrature_point_);
  this->test_set_.SetReflection(this->quadrature_point_,
                                            this->third_quadrature_point_);
  // Second point now has no reflection
  EXPECT_EQ(nullptr, this->test_set_.GetReflection(this->second_quadrature_point_));
  // Third and first now reflect each other
  EXPECT_EQ(this->test_set_.GetReflection(this->quadrature_point_),
            this->third_quadrature_point_);
  EXPECT_EQ(this->test_set_.GetReflection(this->third_quadrature_point_),
            this->quadrature_point_);
}

TYPED_TEST(QuadratureSetTest, SetNewReflectionSecondPoint) {
  this->test_set_.AddPoint(this->third_quadrature_point_);

  this->test_set_.SetReflection(this->quadrature_point_,
                                this->second_quadrature_point_);
  this->test_set_.SetReflection(this->third_quadrature_point_,
                                this->quadrature_point_);
  // Second point now has no reflection
  EXPECT_EQ(nullptr, this->test_set_.GetReflection(this->second_quadrature_point_));
  // Third and first now reflect each other
  EXPECT_EQ(this->test_set_.GetReflection(this->quadrature_point_),
            this->third_quadrature_point_);
  EXPECT_EQ(this->test_set_.GetReflection(this->third_quadrature_point_),
            this->quadrature_point_);
}





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