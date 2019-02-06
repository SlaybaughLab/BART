#ifndef BART_SRC_DATA_FINITE_ELEMENT_GAUSSIAN_H_
#define BART_SRC_DATA_FINITE_ELEMENT_GAUSSIAN_H_

#include <memory>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include "../problem/parameter_types.h"

namespace bart {

namespace data {

template <int dim>
class FiniteElementGaussian {
 public:
  using DiscretizationType = problem::DiscretizationType;
  FiniteElementGaussian(DiscretizationType discretization,
                int polynomial_degree);
  ~FiniteElementGaussian() = default;

  int polynomial_degree() { return polynomial_degree_; };
  dealii::FiniteElement<dim, dim> *finite_element() {
    return finite_element_.get(); };

  int dofs_per_cell() const { return finite_element_->dofs_per_cell; }

  int n_cell_quad_pts() const { return cell_quadrature_->size(); }

  int n_face_quad_pts() const { return face_quadrature_->size(); }

  dealii::FEValues<dim> *finite_element_values() {
    return finite_element_values_.get(); };

  dealii::FEFaceValues<dim> *finite_element_face_values() {
    return finite_element_face_values_.get(); };

  dealii::FEFaceValues<dim> *finite_element_neighbor_face_values() {
    return finite_element_neighbor_face_values_.get(); };
  
  dealii::QGauss<dim> *cell_quadrature() {
    return cell_quadrature_.get(); };

  dealii::QGauss<dim - 1> *face_quadrature() {
    return face_quadrature_.get(); };

 private:
  const int polynomial_degree_;
  std::shared_ptr<dealii::FiniteElement<dim, dim>> finite_element_;
  std::shared_ptr<dealii::FEValues<dim>> finite_element_values_;
  std::shared_ptr<dealii::FEFaceValues<dim>> finite_element_face_values_;
  std::shared_ptr<dealii::FEFaceValues<dim>>
  finite_element_neighbor_face_values_;
  std::shared_ptr<dealii::QGauss<dim>> cell_quadrature_;
  std::shared_ptr<dealii::QGauss<dim - 1>> face_quadrature_;

  std::shared_ptr<dealii::FiniteElement<dim, dim>>
  GetFiniteElement(DiscretizationType discretization);
};

} // namespace data

} // namespace bart 

#endif // BART_SRC_DATA_FINITE_ELEMENT_GAUSSIAN_H_
