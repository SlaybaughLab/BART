#ifndef BART_DOMAIN_FINITE_ELEMENT_HPP_
#define BART_DOMAIN_FINITE_ELEMENT_HPP_

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include "domain/finite_element/finite_element_i.hpp"

namespace bart::domain::finite_element {

/*! \brief Base implementation of the FiniteElement class using Dealii.
 *
 * This class is almost entirley a wrapper for various Dealii objects and many calls are just passed through to the
 * underlying classes.
 *
 * @tparam dim spatial dimension
 */
template<int dim>
class FiniteElement : public FiniteElementI<dim> {
 public:
  using typename FiniteElementI<dim>::DealiiVector, typename FiniteElementI<dim>::Tensor;
  virtual ~FiniteElement() = default;

  // Basic Finite Element data
  [[nodiscard]] auto dofs_per_cell() const -> int override { return finite_element_->dofs_per_cell; }
  [[nodiscard]] auto n_cell_quad_pts() const -> int override { return cell_quadrature_->size(); }
  [[nodiscard]] auto n_face_quad_pts() const -> int override { return face_quadrature_->size(); }

  auto SetCell(const domain::CellPtr<dim> &to_set) -> bool override;
  auto SetFace(const domain::CellPtr<dim> &to_set, const domain::FaceIndex face) -> bool override;

  [[nodiscard]] auto ShapeValue(const int cell_degree_of_freedom,
                                const int cell_quadrature_point) const -> double override {
    return values_->shape_value(cell_degree_of_freedom, cell_quadrature_point);
  }

  [[nodiscard]] auto FaceShapeValue(const int cell_degree_of_freedom,
                                    const int face_quadrature_point) const -> double override {
    return face_values_->shape_value(cell_degree_of_freedom, face_quadrature_point);
  }

  [[nodiscard]] auto ShapeGradient(const int cell_degree_of_freedom,
                                   const int cell_quadrature_point) const -> Tensor override {
    return values_->shape_grad(cell_degree_of_freedom, cell_quadrature_point);
  }

  [[nodiscard]] auto Jacobian(const int cell_quadrature_point) const -> double override {
    return values_->JxW(cell_quadrature_point);
  }

  [[nodiscard]] auto FaceJacobian(const int face_quadrature_point) const -> double override {
    return face_values_->JxW(face_quadrature_point);
  }

  [[nodiscard]] auto FaceNormal() const -> Tensor override {
    return face_values_->normal_vector(0);
  }

  [[nodiscard]] auto ValueAtQuadrature(const DealiiVector& values_at_dofs) const -> std::vector<double> override;
  [[nodiscard]] auto ValueAtFaceQuadrature(const DealiiVector& values_at_dofs) const -> std::vector<double> override;

  // Dependencies
  auto finite_element() -> dealii::FiniteElement<dim, dim>* override { return finite_element_.get(); };
  auto values() -> dealii::FEValues<dim>* override { return values_.get(); };
  auto face_values() -> dealii::FEFaceValues<dim>* override {return face_values_.get(); };
  auto neighbor_face_values() -> dealii::FEFaceValues<dim>* override {return neighbor_face_values_.get(); };
  auto cell_quadrature() -> dealii::QGauss<dim>* { return cell_quadrature_.get(); };
  auto face_quadrature() -> dealii::QGauss<dim - 1>* { return face_quadrature_.get(); };

 protected:
  std::shared_ptr<dealii::FiniteElement<dim, dim>> finite_element_;
  std::shared_ptr<dealii::FEValues<dim>> values_;
  std::shared_ptr<dealii::FEFaceValues<dim>> face_values_;
  std::shared_ptr<dealii::FEFaceValues<dim>> neighbor_face_values_;
  std::shared_ptr<dealii::QGauss<dim>> cell_quadrature_;
  std::shared_ptr<dealii::QGauss<dim - 1>> face_quadrature_;

  bool values_reinit_called_{ false };
  bool face_values_reinit_called_{ false };
  bool neighbor_face_values_reinit_called_{ false };
  using FiniteElementI<dim>::values;
  using FiniteElementI<dim>::face_values;
};

} // namespace bart::domain::finite_element

#endif // BART_DOMAIN_FINITE_ELEMENT_HPP_