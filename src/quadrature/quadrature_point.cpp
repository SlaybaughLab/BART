#include "quadrature/quadrature_point.hpp"

namespace bart::quadrature {

template <int dim>
QuadraturePoint<dim>::QuadraturePoint(std::shared_ptr<quadrature::OrdinateI<dim>> ordinate, const Weight weight)
    : ordinate_(std::move(ordinate)), weight_(weight.get()) {}

template<int dim>
auto QuadraturePoint<dim>::SetOrdinate(const std::shared_ptr<OrdinateI<dim>> ordinate_ptr) -> QuadraturePoint<dim>& {
  ordinate_ = ordinate_ptr;
  return *this;
}

template<int dim>
auto QuadraturePoint<dim>::SetWeight(const Weight weight) -> QuadraturePoint<dim>& {
  weight_ = weight.get();
  return *this;
}

template<int dim>
auto QuadraturePoint<dim>::SetTo(const std::shared_ptr<OrdinateI<dim>> ordinate_ptr,
                                 const Weight weight) -> QuadraturePoint<dim>& {
  SetWeight(weight);
  SetOrdinate(ordinate_ptr);
  return *this;
}

template class QuadraturePoint<1>;
template class QuadraturePoint<2>;
template class QuadraturePoint<3>;

} // namespace bart::quadrature