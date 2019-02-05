#include "finite_element.h"

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_update_flags.h>

namespace bart {

namespace data {

template <int dim>
FiniteElement<dim>::FiniteElement(DiscretizationType discretization,
                                  int polynomial_degree)
    : polynomial_degree_(polynomial_degree) {

  const dealii::UpdateFlags update_flags =
      dealii::update_values | dealii::update_gradients |
      dealii::update_quadrature_points |
      dealii::update_JxW_values;
  
  finite_element_ = GetFiniteElement(discretization);
  
  cell_quadrature_ =
      std::make_shared<dealii::QGauss<dim>>(polynomial_degree + 1);
  
  face_quadrature_ =
      std::make_shared<dealii::QGauss<dim - 1>>(polynomial_degree + 1);

  finite_element_values_ =
      std::make_shared<dealii::FEValues<dim>>(*finite_element_,
                                              *cell_quadrature_,
                                              update_flags);
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
