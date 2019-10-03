#ifndef BART_SRC_QUADRATURE_QUADRATURE_POINT_H_
#define BART_SRC_QUADRATURE_QUADRATURE_POINT_H_

#include "quadrature/quadrature_point_i.h"

namespace bart {

namespace quadrature {

template <int dim>
class QuadraturePoint : public QuadraturePointI<dim> {

};

} // namespace quadrature

} //namespace bart

#endif //BART_SRC_QUADRATURE_QUADRATURE_POINT_H_
