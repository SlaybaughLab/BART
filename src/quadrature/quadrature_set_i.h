#ifndef BART_SRC_QUADRATURE_QUADRATURE_SET_I_H_
#define BART_SRC_QUADRATURE_QUADRATURE_SET_I_H_

#include "quadrature/quadrature_point_i.h"

namespace bart {

namespace quadrature {

/*! \brief Interface for quadrature sets.
 *
 * A quadrature set is a wrapper for a mapping of quadrature points
 * to an index that uniquely defines each of them. Most memeber functions are
 * for accessing points by index or vice versa. The base interface allows for
 * customization of the types of quadrature points.
 *
 * @tparam dim spatial dimension
 */
template <int dim>
class QuadratureSetI {
 public:
  virtual ~QuadratureSetI() = default;
};

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_QUADRATURE_SET_I_H_
