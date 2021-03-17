#ifndef BART_SRC_QUADRATURE_TESTS_ORDINATE_MOCK_HPP_
#define BART_SRC_QUADRATURE_TESTS_ORDINATE_MOCK_HPP_

#include "quadrature/ordinate_i.hpp"

#include "test_helpers/gmock_wrapper.h"

namespace bart::quadrature {

template <int dim>
class OrdinateMock : public OrdinateI<dim> {
 public:
  MOCK_METHOD(OrdinateMock&, set_cartesian_position, (const CartesianPosition<dim>), (override));
  MOCK_METHOD((std::array<double, dim>), cartesian_position, (), (const, override));
  MOCK_METHOD((dealii::Tensor<1, dim>), cartesian_position_tensor, (), (const, override));
  MOCK_METHOD(bool, EquivalenceOperator, (const OrdinateI<dim>&), (const));
  MOCK_METHOD(bool, EquivalenceOperator, ((const std::array<double, dim>)), (const));

  auto operator==(const OrdinateI<dim>& rhs) const -> bool { return EquivalenceOperator(rhs); }
  auto operator!=(const OrdinateI<dim>& rhs) const -> bool { return !EquivalenceOperator(rhs); };
  auto operator==(const std::array<double, dim> rhs) const -> bool { return EquivalenceOperator(rhs); }
  auto operator!=(const std::array<double, dim> rhs) const -> bool { return !EquivalenceOperator(rhs); }
};

} // namespace bart::quadrature

#endif //BART_SRC_QUADRATURE_TESTS_ORDINATE_MOCK_HPP_
