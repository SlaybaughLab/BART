#ifndef BART_SRC_DATA_FINITE_ELEMENT_H_
#define BART_SRC_DATA_FINITE_ELEMENT_H_

#include <memory>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include "../problem/parameter_types.h"

namespace bart {

namespace data {

template <int dim>
class FiniteElement {
 public:
  using DiscretizationType = problem::DiscretizationType;
  FiniteElement(DiscretizationType discretization,
                int polynomial_degree);
  ~FiniteElement() = default;

  int polynomial_degree() { return polynomial_degree_; };
  dealii::FiniteElement<dim, dim> *finite_element() {
    return finite_element_.get(); };

  dealii::FEValues<dim> *finite_element_values_() {
    return finite_element_values_.get(); };

  dealii::QGauss<dim> *cell_quadrature() {
    return cell_quadrature_.get(); };

  dealii::QGauss<dim - 1> *face_quadrature() {
    return face_quadrature_.get(); };

 private:
  const int polynomial_degree_;
  std::shared_ptr<dealii::FiniteElement<dim, dim>> finite_element_;
  std::shared_ptr<dealii::FEValues<dim>> finite_element_values_;
  std::shared_ptr<dealii::QGauss<dim>> cell_quadrature_;
  std::shared_ptr<dealii::QGauss<dim - 1>> face_quadrature_;


  std::shared_ptr<dealii::FiniteElement<dim, dim>>
  GetFiniteElement(DiscretizationType discretization);
};

} // namespace data

} // namespace bart 

#endif // BART_SRC_DATA_FINITE_ELEMENT_H_
