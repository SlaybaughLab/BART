#include "quadrature/quadrature_point.h"

#include "quadrature/tests/ordinate_mock.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

/* Tests to verify the operation of the default implementation of the
 * QuadraturePoint class. Mocking of an ordinate object is required.
 */
template <typename DimensionWrapper>
class QuadraturePointTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  std::shared_ptr<quadrature::OrdinateI<dim>> ordinate_ptr;
  void SetUp() override;
};

// SetUp: instantiates mock ordinate
template <typename DimensionWrapper>
void QuadraturePointTest<DimensionWrapper>::SetUp() {
  ordinate_ptr = std::make_shared<quadrature::OrdinateMock<dim>>();
}

TYPED_TEST_CASE(QuadraturePointTest, bart::testing::AllDimensions);

/* Default constructor should initialize pointer as null and weight as 0.
 * Constructor taking ordinate and weight should set them properly. Also tests
 * getters to verify.
 */
TYPED_TEST(QuadraturePointTest, ConstructorAndGetters) {
  constexpr int dim = this->dim;
  const double weight = 1.45;

  quadrature::QuadraturePoint<dim> default_point;
  EXPECT_EQ(default_point.ordinate(), nullptr);
  EXPECT_EQ(default_point.weight(), 0);

  quadrature::QuadraturePoint<dim> new_point(this->ordinate_ptr,
                                             quadrature::Weight(weight));

  EXPECT_EQ(this->ordinate_ptr.use_count(), 2);
  EXPECT_EQ(new_point.ordinate(), this->ordinate_ptr);
  EXPECT_EQ(new_point.weight(), weight);
}

// Setters should set quadrature point ordinate and weight properly
TYPED_TEST(QuadraturePointTest, Setters) {
  constexpr int dim = this->dim;
  const double weight = 1.45, second_weight = 2.1;

  quadrature::QuadraturePoint<dim> test_point;
  auto second_point_ptr = std::make_shared<quadrature::OrdinateMock<dim>>();

  test_point.SetOrdinate(this->ordinate_ptr);
  test_point.SetWeight(quadrature::Weight(weight));

  EXPECT_GT(this->ordinate_ptr.use_count(), 1);
  EXPECT_EQ(test_point.ordinate(), this->ordinate_ptr);
  EXPECT_EQ(test_point.weight(), weight);

  test_point.SetOrdinate(second_point_ptr);
  test_point.SetWeight(quadrature::Weight(second_weight));

  EXPECT_GT(second_point_ptr.use_count(), 1);
  EXPECT_EQ(test_point.ordinate(), second_point_ptr);
  EXPECT_EQ(test_point.weight(), second_weight);

  test_point.SetTo(this->ordinate_ptr, quadrature::Weight(weight));
  EXPECT_GT(this->ordinate_ptr.use_count(), 1);
  EXPECT_EQ(test_point.ordinate(), this->ordinate_ptr);
  EXPECT_EQ(test_point.weight(), weight);
}

// Cartesian position should properly return the ordinate position.
TYPED_TEST(QuadraturePointTest, CartesianPositionAndTensor) {
  constexpr int dim = this->dim;

  quadrature::QuadraturePoint<dim> test_point;

  test_point.SetOrdinate(this->ordinate_ptr);

  std::array<double, dim> position;
  position.fill(1);

  auto mock_ptr = dynamic_cast<quadrature::OrdinateMock<dim>*>(this->ordinate_ptr.get());

  EXPECT_CALL(*mock_ptr, cartesian_position())
      .WillOnce(::testing::Return(position));

  dealii::Tensor<1, dim> tensor_position;
  for (int i = 0; i < dim; ++i)
    tensor_position[i] = 1;

  EXPECT_CALL(*mock_ptr, cartesian_position_tensor())
      .WillOnce(::testing::Return(tensor_position));

  EXPECT_EQ(position, test_point.cartesian_position());
  auto tensor = test_point.cartesian_position_tensor();

  for (int i = 0; i < dim; ++i)
    EXPECT_EQ(tensor[i], 1);
}

} // namespace