#ifndef BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SET_I_H_
#define BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SET_I_H_

#include "quadrature/angular/angular_quadrature_types.h"

namespace bart {

namespace quadrature {

namespace angular {

/*! \brief Provides the interface for angular quadrature sets.
 *
 * Angular colocation methods provide sets of
 *
 */
class AngularQuadratureSetI {
 public:
  virtual ~AngularQuadratureSetI() = default;

};

} // namespace angular

} // namespace quadrature

} // namespace bart

#endif // BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SET_I_H_