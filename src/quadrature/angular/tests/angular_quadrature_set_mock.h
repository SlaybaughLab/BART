#ifndef BART_SRC_QUADRATURE_ANGULAR_TESTS_ANGULAR_QUADRATURE_SET_MOCK_H_
#define BART_SRC_QUADRATURE_ANGULAR_TESTS_ANGULAR_QUADRATURE_SET_MOCK_H_

#include "quadrature/angular/angular_quadrature_set_i.h"

#include "quadrature/angular/angular_quadrature_types.h"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace quadrature {

namespace angular {

template <int dim>
class AngularQuadratureSetMock : public AngularQuadratureSetI<dim> {
 public:
  using typename AngularQuadratureSetI<dim>::AngleIndex;
  MOCK_CONST_METHOD0_T(quadrature_points_map,
                       std::map<AngleIndex, QuadraturePoint<dim>>());
  MOCK_CONST_METHOD0_T(quadrature_points,
                       std::vector<QuadraturePoint<dim>>());
  MOCK_CONST_METHOD0_T(quadrature_weights, std::vector<Weight>());
  MOCK_CONST_METHOD0_T(total_quadrature_points, int());
};

} // namespace angular

} // namespace quadrature

} // namespace bart

#endif // BART_SRC_QUADRATURE_ANGULAR_TESTS_ANGULAR_QUADRATURE_SET_MOCK_H_