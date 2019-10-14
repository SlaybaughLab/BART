#include "quadrature/quadrature_point.h"

namespace bart {

namespace quadrature {

template <int dim>
QuadraturePoint<dim>::QuadraturePoint(
    std::shared_ptr<quadrature::OrdinateI<dim>> ordinate, Weight weight)
    : ordinate_(ordinate),
      weight_(weight.get()) {}

template<int dim>
QuadraturePoint<dim> &QuadraturePoint<dim>::SetOrdinate(
    const std::shared_ptr<OrdinateI<dim>>& ordinate_ptr) {
  ordinate_ = ordinate_ptr;
  return *this;
}

template<int dim>
QuadraturePoint<dim> &QuadraturePoint<dim>::SetWeight(const Weight weight) {
  weight_ = weight.get();
  return *this;
}
template<int dim>
QuadraturePoint<dim> &QuadraturePoint<dim>::SetTo(
    const std::shared_ptr<OrdinateI<dim>> &ordinate_ptr, const Weight weight) {
  SetWeight(weight);
  SetOrdinate(ordinate_ptr);
  return *this;
}

template class QuadraturePoint<1>;
template class QuadraturePoint<2>;
template class QuadraturePoint<3>;

} // namespace quadrature

} // namespace bart
