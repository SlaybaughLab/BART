#include "finite_element.h"

namespace bart {

namespace data {

template <int dim>
FiniteElement<dim>::FiniteElement(bart::problem::DiscretizationType discretization,
                                  int polynomial_degree)
    : polynomial_degree_(polynomial_degree) {
  // Generate appropriate discretization
  
}

template class FiniteElement<1>;
template class FiniteElement<2>;
template class FiniteElement<3>;

} // namespace data

} // namespace bart 
