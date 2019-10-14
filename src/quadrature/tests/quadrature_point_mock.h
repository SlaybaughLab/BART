#ifndef BART_SRC_QUADRATURE_TESTS_QUADRATURE_POINT_MOCK_H_
#define BART_SRC_QUADRATURE_TESTS_QUADRATURE_POINT_MOCK_H_

#include "quadrature/quadrature_point_i.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace quadrature {

template <int dim>
class QuadraturePointMock : public QuadraturePointI<dim> {
 public:
  MOCK_METHOD2_T(SetTo, QuadraturePointMock&(const std::shared_ptr<OrdinateI<dim>>&,
      const quadrature::Weight));
  MOCK_METHOD1_T(SetOrdinate, QuadraturePointMock&(const std::shared_ptr<OrdinateI<dim>>&));
  MOCK_METHOD1_T(SetWeight, QuadraturePointMock&(const quadrature::Weight));
  MOCK_CONST_METHOD0_T(ordinate, std::shared_ptr<OrdinateI<dim>>());
  MOCK_CONST_METHOD0_T(weight, double());
};

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_TESTS_QUADRATURE_POINT_MOCK_H_
