#include "domain/finite_element/finite_element_gaussian.hpp"

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_update_flags.h>

namespace bart::domain::finite_element {

template <int dim>
FiniteElementGaussian<dim>::FiniteElementGaussian(DiscretizationType discretization, int polynomial_degree)
    : polynomial_degree_(polynomial_degree) {
  std::string description{"deal.II Gaussian, " + std::to_string(dim) + "D, "};

  const auto update_flags = dealii::update_values | dealii::update_gradients | dealii::update_quadrature_points |
                            dealii::update_JxW_values;
  const auto face_update_flags = dealii::update_values | dealii::update_gradients | dealii::update_quadrature_points |
      dealii::update_JxW_values | dealii::update_normal_vectors;
  
  finite_element_ = GetFiniteElement(discretization);
  cell_quadrature_ = std::make_shared<dealii::QGauss<dim>>(polynomial_degree + 1);
  face_quadrature_ = std::make_shared<dealii::QGauss<dim - 1>>(polynomial_degree + 1);
  values_ = std::make_shared<dealii::FEValues<dim>>(*finite_element_, *cell_quadrature_, update_flags);
  face_values_ = std::make_shared<dealii::FEFaceValues<dim>>(*finite_element_, *face_quadrature_, face_update_flags);
  
  if (discretization == DiscretizationType::kDiscontinuousFEM) {
    neighbor_face_values_ = 
        std::make_shared<dealii::FEFaceValues<dim>>(*finite_element_, *face_quadrature_, face_update_flags);
    description += "Discontinuous, ";
  } else {
    description += "Continuous, ";
  }
  description += "Q = " + std::to_string(polynomial_degree);
  this->set_description(description, utility::DefaultImplementation(true));
}

template <int dim>
auto FiniteElementGaussian<dim>::GetFiniteElement(DiscretizationType discretization)
-> std::shared_ptr<dealii::FiniteElement<dim, dim>> {
  switch (discretization) {
    case DiscretizationType::kContinuousFEM: {
      return std::make_shared<dealii::FE_Q<dim>>(polynomial_degree_);
    }
    case DiscretizationType::kDiscontinuousFEM: {
      return std::make_shared<dealii::FE_DGQ<dim>>(polynomial_degree_);
    }
    default: {
      AssertThrow(false, dealii::ExcMessage("Cannot build FiniteElementGaussian object with discretization type None"));
      break;
    }
  }
}

template class FiniteElementGaussian<1>;
template class FiniteElementGaussian<2>;
template class FiniteElementGaussian<3>;

} // namespace bart::domain::finite_element
