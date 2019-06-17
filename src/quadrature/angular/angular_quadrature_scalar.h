#ifndef BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SCALAR_H_
#define BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SCALAR_H_

#include "quadrature/angular/angular_quadrature_set.h"

namespace bart {

namespace quadrature {

namespace angular {

/*! \brief Provides the angular quadrature for scalar problems.
 *
 * For scalar problems, there is no angular dependence, but many functions
 * require an angular quadrature for generalization. This class provides one
 * possible "scalar angular quadrature":
 *
 * \f[
 * (w_0, \hat{\Omega}_0) = (1, \vec{0})\;,
 * \f]
 * where \f$\vec{0} \in \mathbb{R}^n\f$ and \f$n\f$ is the angular dimension of
 * the problem.
 *
 */
template <int dim>
class AngularQuadratureScalar : public AngularQuadratureSet<dim> {
 public:
  AngularQuadratureScalar() = default;
  virtual ~AngularQuadratureScalar() = default;

};

} // namespace angular

} // namespace quadrature

} // namespace bart

#endif // BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SCALAR_H_