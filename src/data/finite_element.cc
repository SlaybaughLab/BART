#include "finite_element.h"

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>

namespace bart {

namespace data {

template <int dim>
FiniteElement<dim>::FiniteElement(DiscretizationType discretization,
                                  int polynomial_degree)
    : polynomial_degree_(polynomial_degree) {
  
  finite_element_ = GetFiniteElement(discretization);
}

template <int dim>
std::shared_ptr<dealii::FiniteElement<dim, dim>>
FiniteElement<dim>::GetFiniteElement(DiscretizationType discretization) {
  switch (discretization) {
    case DiscretizationType::kContinuousFEM: {
      return std::make_shared<dealii::FE_Q<dim>>(polynomial_degree_);
      break;
    }
    case DiscretizationType::kDiscontinuousFEM:{
      return std::make_shared<dealii::FE_DGQ<dim>>(polynomial_degree_);
      break;
    }
    default: {
      AssertThrow(false,
                  dealii::ExcMessage("Cannot build FiniteElement object with discretization type None"));
      break;
    }
  }
}

template class FiniteElement<1>;
template class FiniteElement<2>;
template class FiniteElement<3>;

} // namespace data

} // namespace bart 
