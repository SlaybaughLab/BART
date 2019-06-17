#ifndef BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SCALAR_H_
#define BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SCALAR_H_

#include "quadrature/angular/angular_quadrature_set_i.h"

namespace bart {

namespace quadrature {

namespace angular {

class AngularQuadratureScalar : public AngularQuadratureSetI {
 public:
  virtual ~AngularQuadratureScalar() = default;
};

} // namespace angular

} // namespace quadrature

} // namespace bart

#endif // BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SCALAR_H_