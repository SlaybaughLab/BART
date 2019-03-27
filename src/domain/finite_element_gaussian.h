#ifndef BART_SRC_DOMAIN_FINITE_ELEMENT_GAUSSIAN_H_
#define BART_SRC_DOMAIN_FINITE_ELEMENT_GAUSSIAN_H_

#include <memory>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

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
  
  dealii::FiniteElement<dim, dim> *finite_element() override {
    return finite_element_.get(); };

  int dofs_per_cell() const override { return finite_element_->dofs_per_cell; }

  int n_cell_quad_pts() const override { return cell_quadrature_->size(); }

  int n_face_quad_pts() const override { return face_quadrature_->size(); }

  dealii::FEValues<dim> *values() override { return values_.get(); };

  dealii::FEFaceValues<dim> *face_values() override {
    return face_values_.get(); };

  dealii::FEFaceValues<dim> *neighbor_face_values() override {
    return neighbor_face_values_.get(); };
  
  dealii::QGauss<dim> *cell_quadrature() { return cell_quadrature_.get(); };

  dealii::QGauss<dim - 1> *face_quadrature() { return face_quadrature_.get(); };

 private:
  const int polynomial_degree_;
  std::shared_ptr<dealii::FiniteElement<dim, dim>> finite_element_;
  std::shared_ptr<dealii::FEValues<dim>> values_;
  std::shared_ptr<dealii::FEFaceValues<dim>> face_values_;
  std::shared_ptr<dealii::FEFaceValues<dim>> neighbor_face_values_;
  std::shared_ptr<dealii::QGauss<dim>> cell_quadrature_;
  std::shared_ptr<dealii::QGauss<dim - 1>> face_quadrature_;

  std::shared_ptr<dealii::FiniteElement<dim, dim>>
  GetFiniteElement(DiscretizationType discretization);
};

} // namespace domain

} // namespace bart 

#endif // BART_SRC_DOMAIN_FINITE_ELEMENT_GAUSSIAN_H_
