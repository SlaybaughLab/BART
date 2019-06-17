#ifndef BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SET_I_H_
#define BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SET_I_H_

#include "quadrature/angular/angular_quadrature_types.h"
#include "data/system/system_types.h"

#include <map>
#include <vector>

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
 * where \f$n\f$ is the angular dimension of the problem, and \f$N\f$ is the
 * total number of quadrature points. In general, the set is chosen such that
 * it accurately integrates specific functions in specific domains.
 *
 * @tparam dim angular dimension \f$n\f$ of the problem
 */
template <int dim>
class AngularQuadratureSetI {
 public:
  using AngleIndex = data::system::AngleIndex;
  virtual ~AngularQuadratureSetI() = default;

  /*! \brief Map of quadrature points to system angle indices. */
  virtual std::map<AngleIndex, QuadraturePoint<dim>> quadrature_points_map() const = 0;

  /*! \brief Returns a vector holding the quadrature points
   *
   * Position in the vector corresponds to the AngleIndex used by solutions.
   *
   * */
  virtual std::vector<QuadraturePoint<dim>> quadrature_points() const = 0;
  /*! \brief Returns a vector holding the quadrature weights
   *
   * Position in the vector corresponds to the AngleIndex used by solutions.
   *
   * */
  virtual std::vector<Weight> quadrature_weights() const = 0;

  virtual int total_quadrature_points() const = 0;

};

} // namespace angular

} // namespace quadrature

} // namespace bart

#endif // BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SET_I_H_