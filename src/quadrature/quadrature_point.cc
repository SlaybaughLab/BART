#include "quadrature/quadrature_point.h"

namespace bart {

namespace quadrature {

template <int dim>
QuadraturePoint<dim>::QuadraturePoint(
    std::shared_ptr<quadrature::OrdinateI<dim>> ordinate, Weight weight)
    : ordinate_(ordinate),
      weight_(weight.get()) {}

template class QuadraturePoint<1>;
template class QuadraturePoint<2>;
template class QuadraturePoint<3>;

} // namespace quadrature

} // namespace bart
