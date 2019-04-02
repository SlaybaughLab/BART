#ifndef BART_SRC_DOMAIN_FINITE_ELEMENT_GAUSSIAN_H_
#define BART_SRC_DOMAIN_FINITE_ELEMENT_GAUSSIAN_H_

#include <memory>

#include "domain/finite_element.h"
#include "problem/parameter_types.h"

namespace bart {

namespace domain {

template <int dim>
class FiniteElementGaussian : public FiniteElement<dim>{
 public:
  using DiscretizationType = problem::DiscretizationType;
  FiniteElementGaussian(DiscretizationType discretization,
                int polynomial_degree);
  ~FiniteElementGaussian() = default;

  int polynomial_degree() const override { return polynomial_degree_; };

 private:
  const int polynomial_degree_;
  using FiniteElement<dim>::finite_element_;
  using FiniteElement<dim>::values_;
  using FiniteElement<dim>::face_values_;
  using FiniteElement<dim>::neighbor_face_values_;
  using FiniteElement<dim>::cell_quadrature_;
  using FiniteElement<dim>::face_quadrature_;

  std::shared_ptr<dealii::FiniteElement<dim, dim>>
  GetFiniteElement(DiscretizationType discretization);
};

} // namespace domain

} // namespace bart 

#endif // BART_SRC_DOMAIN_FINITE_ELEMENT_GAUSSIAN_H_
