#ifndef BART_SRC_QUADRATURE_ANGULAR_TESTS_ANGULAR_QUADRATURE_GENERATOR_MOCK_H_
#define BART_SRC_QUADRATURE_ANGULAR_TESTS_ANGULAR_QUADRATURE_GENERATOR_MOCK_H_

#include "quadrature/quadrature_generator_i.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace quadrature {

template <int dim>
class QuadratureGeneratorMock : QuadratureGeneratorI<dim> {
 public:
  MOCK_CONST_METHOD0_T(GenerateSet, std::vector<std::pair<CartesianPosition<dim>, Weight>>());
  MOCK_CONST_METHOD0_T(order, int());
};

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_ANGULAR_TESTS_ANGULAR_QUADRATURE_GENERATOR_MOCK_H_
