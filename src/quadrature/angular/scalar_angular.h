#ifndef BART_SRC_QUADRATURE_ANGULAR_SCALAR_ANGULAR_H_
#define BART_SRC_QUADRATURE_ANGULAR_SCALAR_ANGULAR_H_

#include "quadrature/quadrature_generator_i.h"

namespace bart {

namespace quadrature {

namespace angular {

template <int dim>
class ScalarAngular : public QuadratureGeneratorI<dim> {
 public:
  ScalarAngular() {};
  std::vector<std::pair<CartesianPosition<dim>, Weight>> GenerateSet() const override {
    std::vector<std::pair<CartesianPosition<dim>, Weight>> generated_set;
    std::array<double, dim> null_position{0};
    generated_set.emplace_back(std::make_pair(CartesianPosition<dim>(null_position),
                                              Weight(1.0)));
    return generated_set;
  };
  int order() const override { return 0; };
};

} // namespace angular

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_ANGULAR_SCALAR_ANGULAR_H_
