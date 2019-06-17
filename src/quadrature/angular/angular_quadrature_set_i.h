#ifndef BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SET_I_H_
#define BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SET_I_H_

#include "quadrature/angular/angular_quadrature_types.h"

namespace bart {

namespace quadrature {

namespace angular {

/*! \brief Provides the interface for angular quadrature sets.
 *
 * Sets provide a quadrature set, defined as a set of positions (ordinates) and
 * weights. The set is defined as:
 *
 * \f[
 * \{(w_i, \hat{\Omega}_i) \mid w \in \mathbb{R}, \hat{\Omega} \in \mathbb{R}^n, 0 \leq i \leq N\}\;,
 * \f]
 *
 * where \f$n\f$ is the spatial dimension of the problem, and \f$N\f$ is the
 * total number of quadrature points. In general, the set is chosen such that
 * it accurately integrates specific functions in specific domains.
 *
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