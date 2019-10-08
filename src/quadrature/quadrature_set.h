#ifndef BART_SRC_QUADRATURE_QUADRATURE_SET_H_
#define BART_SRC_QUADRATURE_QUADRATURE_SET_H_

#include "quadrature/quadrature_set_i.h"

namespace bart {

namespace quadrature {

/*! \brief Default implementation of the QuadratureSet.
 *
 * @tparam dim spatial dimension
 */
template <int dim>
class QuadratureSet : public QuadratureSetI<dim> {
 public:
  virtual ~QuadratureSet() = default;

};

} // namespace quadrature

} //namespace bart

#endif //BART_SRC_QUADRATURE_QUADRATURE_SET_H_
