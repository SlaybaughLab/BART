#include "domain/finite_element/finite_element_factory.h"

#include "domain/finite_element/finite_element_gaussian.h"

namespace bart {

namespace domain {

namespace finite_element {

template<int dim>
std::unique_ptr<FiniteElementI<dim>>
    FiniteElementFactory<dim>::MakeFiniteElement(
        const problem::DiscretizationType discretization_type,
        const int polynomial_degree,
        const FiniteElementImpl finite_element_type) {
  std::unique_ptr<FiniteElementI<dim>> return_ptr = nullptr;

  if (finite_element_type == FiniteElementImpl::kGaussian) {
    return_ptr = std::move(
        std::make_unique<FiniteElementGaussian<dim>>(
            discretization_type, polynomial_degree));
  }

  return return_ptr;
}

template class FiniteElementFactory<1>;
template class FiniteElementFactory<2>;
template class FiniteElementFactory<3>;

} // namespace finite_element

} // namespace domain

} //namespace bart
