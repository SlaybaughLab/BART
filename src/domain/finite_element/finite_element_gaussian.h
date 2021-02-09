#ifndef BART_SRC_DOMAIN_FINITE_ELEMENT_GAUSSIAN_H_
#define BART_SRC_DOMAIN_FINITE_ELEMENT_GAUSSIAN_H_

#include <memory>

#include "domain/finite_element/finite_element.hpp"
#include "problem/parameter_types.hpp"

namespace bart {

namespace domain {

namespace finite_element {

/*! \brief Provides finite element information using a Gaussian cell quadrature.
 *
 * For a continuous FEM basis, this object uses the FE_Q<dim> object, and
 * FE_DGQ for a discontinuous basis.
 *
 *
 * \tparam dim dimension.
 */

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

} // namespace finite_element

} // namespace domain

} // namespace bart 

#endif // BART_SRC_DOMAIN_FINITE_ELEMENT_GAUSSIAN_H_
