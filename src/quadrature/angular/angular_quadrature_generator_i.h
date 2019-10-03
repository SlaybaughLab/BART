#ifndef BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_GENERATOR_I_H_
#define BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_GENERATOR_I_H_

#include "quadrature/quadrature_types.h"

#include <utility>
#include <vector>

namespace bart {

namespace quadrature {

namespace angular {

template <int dim>
class AngularQuadratureGeneratorI {
 public:
  virtual ~AngularQuadratureGeneratorI() = default;
  virtual std::vector<std::pair<CartesianPosition<dim>, Weight>> GenerateSet() const = 0;
  virtual int order() const = 0;

};

} // namespace angular

} // namespace quadrature

} // namespace bart




#endif //BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_GENERATOR_I_H_
