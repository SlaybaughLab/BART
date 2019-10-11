#include "quadrature/quadrature_set.h"

#include <memory>

#include "quadrature/tests/quadrature_point_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class QuadratureSetTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;

  quadrature::QuadratureSet<dim> test_set_;
  std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point_,
      second_quadrature_point_;

  void SetUp() override {
    quadrature_point_ =
        std::make_shared<quadrature::QuadraturePointMock<dim>>();
    second_quadrature_point_ =
        std::make_shared<quadrature::QuadraturePointMock<dim>>();
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