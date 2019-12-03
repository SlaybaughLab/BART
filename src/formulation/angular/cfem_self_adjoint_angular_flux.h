#ifndef BART_SRC_FORMULATION_ANGULAR_CFEM_SELF_ADJOINT_ANGULAR_FLUX_H_
#define BART_SRC_FORMULATION_ANGULAR_CFEM_SELF_ADJOINT_ANGULAR_FLUX_H_

#include "data/cross_sections.h"
#include "domain/finite_element/finite_element_i.h"
#include "formulation/angular/cfem_self_adjoint_angular_flux_i.h"
#include "quadrature/quadrature_set_i.h"

#include <memory>

namespace bart {

namespace formulation {

namespace angular {

template <int dim>
class CFEMSelfAdjointAngularFlux : public CFEMSelfAdjointAngularFluxI<dim> {
 public:
  using typename CFEMSelfAdjointAngularFluxI<dim>::InitializationToken;

  CFEMSelfAdjointAngularFlux(
      std::shared_ptr<domain::finite_element::FiniteElementI<dim>>,
      std::shared_ptr<data::CrossSections>,
      std::shared_ptr<quadrature::QuadratureSetI<dim>>);

  InitializationToken Initialize(const formulation::CellPtr<dim>&) override;

  std::vector<double> OmegaDotGradient(
      int cell_quadrature_point,
      quadrature::QuadraturePointIndex angular_index) const {
    return omega_dot_gradient_.at({cell_quadrature_point, angular_index.get()});
  }

  FullMatrix OmegaDotGradientSquared(
      int cell_quadrature_point,
      quadrature::QuadraturePointIndex angular_index) const {
    return omega_dot_gradient_squared_.at(cell_quadrature_point)
        .at(angular_index.get()); }

  // Dependency getters
  domain::finite_element::FiniteElementI<dim>* finite_element_ptr() const {
    return finite_element_ptr_.get(); }
  data::CrossSections* cross_sections_ptr() const {
    return cross_sections_ptr_.get(); }
  quadrature::QuadratureSetI<dim>* quadrature_set_ptr() const {
    return quadrature_set_ptr_.get(); }

  std::map<int, FullMatrix> shape_squared() const {
    return shape_squared_; }

 protected:
  // Dependencies
  std::shared_ptr<domain::finite_element::FiniteElementI<dim>> finite_element_ptr_;
  std::shared_ptr<data::CrossSections> cross_sections_ptr_;
  std::shared_ptr<quadrature::QuadratureSetI<dim>> quadrature_set_ptr_;
  // Geometric properties
  const int cell_degrees_of_freedom_ = 0; //!< Degrees of freedom per cell
  const int cell_quadrature_points_ = 0; //!< Quadrature points per cell
  const int face_quadrature_points_ = 0; //!< Quadrature points per face
  // Precalculated matrices and vectors
  using CellQuadratureIndex = int;
  using AngleIndex = int;
  std::map<std::pair<CellQuadratureIndex, AngleIndex>,
           std::vector<double>> omega_dot_gradient_;
  std::map<int, std::map<int, FullMatrix>> omega_dot_gradient_squared_;
  std::map<CellQuadratureIndex, FullMatrix> shape_squared_ = {};

};

} // namespace angular

} // namespace formulation

} //namespace bart

#endif //BART_SRC_FORMULATION_ANGULAR_CFEM_SELF_ADJOINT_ANGULAR_FLUX_H_
