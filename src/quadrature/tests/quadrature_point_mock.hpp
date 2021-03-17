#ifndef BART_SRC_QUADRATURE_TESTS_QUADRATURE_POINT_MOCK_HPP_
#define BART_SRC_QUADRATURE_TESTS_QUADRATURE_POINT_MOCK_HPP_

#include "quadrature/quadrature_point_i.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart::quadrature {

template <int dim>
class QuadraturePointMock : public QuadraturePointI<dim> {
 public:
  MOCK_METHOD(QuadraturePointMock&, SetTo, (const std::shared_ptr<OrdinateI<dim>>&, const quadrature::Weight), (override));
  MOCK_METHOD(QuadraturePointMock&, SetOrdinate, (const std::shared_ptr<OrdinateI<dim>>&), (override));
  MOCK_METHOD(QuadraturePointMock&, SetWeight, (const quadrature::Weight), (override));
  MOCK_METHOD(std::shared_ptr<OrdinateI<dim>>, ordinate, (), (const, override));
  MOCK_METHOD(double, weight, (), (const, override));
  MOCK_METHOD((std::array<double, dim>), cartesian_position, (), (const, override));
  MOCK_METHOD((dealii::Tensor<1, dim>), cartesian_position_tensor, (), (const, override));
};

} // namespace bart::quadrature

#endif //BART_SRC_QUADRATURE_TESTS_QUADRATURE_POINT_MOCK_HPP_
