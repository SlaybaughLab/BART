#ifndef BART_SRC_QUADRATURE_TESTS_ORDINATE_MOCK_H_
#define BART_SRC_QUADRATURE_TESTS_ORDINATE_MOCK_H_

#include "quadrature/ordinate_i.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace quadrature {

template <int dim>
class OrdinateMock : public OrdinateI<dim> {
 public:
  MOCK_METHOD1_T(set_cartesian_position, OrdinateMock&(const CartesianPosition<dim>));
  MOCK_CONST_METHOD0_T(cartesian_position, std::array<double, dim>());
  MOCK_CONST_METHOD0_T(cartesian_position_tensor, dealii::Tensor<1, dim>());

  MOCK_CONST_METHOD1_T(EquivalenceOperator, bool(const OrdinateI<dim>&));
  MOCK_CONST_METHOD1_T(EquivalenceOperator, bool(const std::array<double, dim>));

  bool operator==(const OrdinateI<dim>& rhs) const {
    return EquivalenceOperator(rhs);
  }
  bool operator!=(const OrdinateI<dim>& rhs) const {
    return !EquivalenceOperator(rhs);
  };
  bool operator==(const std::array<double, dim> rhs) const {
    return EquivalenceOperator(rhs);
  }
  bool operator!=(const std::array<double, dim> rhs) const {
    return !EquivalenceOperator(rhs);
  }
};

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_TESTS_ORDINATE_MOCK_H_
